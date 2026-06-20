// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "lumpedportoperator.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <fmt/ranges.h>
#include "fem/coefficient.hpp"
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "linalg/operator.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

namespace palace
{

using namespace std::complex_literals;

namespace
{

using Point = LumpedPortData::Point;
using TerminalEdge = LumpedPortData::TerminalEdge;

Point GetVertexPoint(const mfem::ParMesh &mesh, int vertex)
{
  const auto *x = mesh.GetVertex(vertex);
  return {x[0], x[1], x[2]};
}

Point Subtract(const Point &a, const Point &b)
{
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

double Dot(const Point &a, const Point &b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double Norm(const Point &a)
{
  return std::sqrt(Dot(a, a));
}

double Distance(const Point &a, const Point &b)
{
  return Norm(Subtract(a, b));
}

std::vector<int> UniqueAttributes(const mfem::Array<int> &attrs)
{
  std::vector<int> unique_attrs;
  unique_attrs.reserve(attrs.Size());
  for (int i = 0; i < attrs.Size(); i++)
  {
    if (std::find(unique_attrs.begin(), unique_attrs.end(), attrs[i]) == unique_attrs.end())
    {
      unique_attrs.push_back(attrs[i]);
    }
  }
  return unique_attrs;
}

void PrintBoundaryAdjacencyDiagnostics(const mfem::ParMesh &mesh,
                                       const mfem::Array<int> &attrs,
                                       const std::string &label)
{
  const auto unique_attrs = UniqueAttributes(attrs);
  if (unique_attrs.empty())
  {
    return;
  }

  const int domain_attr_max = mesh.attributes.Size() ? mesh.attributes.Max() : 0;
  const int stride = domain_attr_max + 1;
  std::vector<long long> counts(unique_attrs.size() * stride, 0);
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    const int bdr_attr = mesh.GetBdrAttribute(be);
    const auto it = std::find(unique_attrs.begin(), unique_attrs.end(), bdr_attr);
    if (it == unique_attrs.end())
    {
      continue;
    }

    int elem_id = -1, info = 0;
    mesh.GetBdrElementAdjacentElement(be, elem_id, info);
    int domain_attr = 0;
    if (elem_id >= 0)
    {
      domain_attr = mesh.GetAttribute(elem_id);
      if (domain_attr < 0 || domain_attr > domain_attr_max)
      {
        domain_attr = 0;
      }
    }

    const std::size_t attr_idx = std::distance(unique_attrs.begin(), it);
    counts[attr_idx * stride + domain_attr]++;
  }
  Mpi::GlobalSum(static_cast<int>(counts.size()), counts.data(), mesh.GetComm());

  fmt::memory_buffer buffer{};
  auto out = fmt::appender{buffer};
  fmt::format_to(out, "\nBoundary adjacency diagnostics for {}:\n", label);
  for (std::size_t attr_idx = 0; attr_idx < unique_attrs.size(); attr_idx++)
  {
    long long total = 0;
    for (int domain_attr = 0; domain_attr <= domain_attr_max; domain_attr++)
    {
      total += counts[attr_idx * stride + domain_attr];
    }
    fmt::format_to(out, " Boundary attr {:d}: total elements = {:d}",
                   unique_attrs[attr_idx], total);
    for (int domain_attr = 0; domain_attr <= domain_attr_max; domain_attr++)
    {
      const long long count = counts[attr_idx * stride + domain_attr];
      if (count == 0)
      {
        continue;
      }
      if (domain_attr == 0)
      {
        fmt::format_to(out, ", adjacent attr <none/invalid> = {:d}", count);
      }
      else
      {
        fmt::format_to(out, ", adjacent attr {:d} = {:d}", domain_attr, count);
      }
    }
    fmt::format_to(out, "\n");
  }
  Mpi::Print("{}", fmt::to_string(buffer));
}

std::pair<double, double> SegmentCoordinateAndDistance(const Point &x,
                                                       const TerminalEdge &edge)
{
  const auto ab = Subtract(edge[1], edge[0]);
  const auto ax = Subtract(x, edge[0]);
  const double length2 = Dot(ab, ab);
  MFEM_VERIFY(length2 > 0.0, "\"TerminalEdges\" entry has coincident endpoints!");
  const double t = Dot(ax, ab) / length2;
  const Point projected{edge[0][0] + t * ab[0], edge[0][1] + t * ab[1],
                        edge[0][2] + t * ab[2]};
  return {t, Distance(x, projected)};
}

LumpedPortData::TerminalEdgeChain FindTerminalEdgeChain(const mfem::ParMesh &mesh,
                                                        const TerminalEdge &edge)
{
  constexpr double distance_tol = 1.0e-8;
  constexpr double coordinate_tol = 1.0e-8;
  struct Candidate
  {
    int edge;
    double t_min;
    double t_max;
    double length;
    double sign;
  };
  std::vector<Candidate> candidates;
  mfem::Array<int> vertices;
  for (int edge_idx = 0; edge_idx < mesh.GetNEdges(); edge_idx++)
  {
    mesh.GetEdgeVertices(edge_idx, vertices);
    MFEM_VERIFY(vertices.Size() == 2, "Expected mesh edge to have two vertices!");
    const auto first = GetVertexPoint(mesh, vertices[0]);
    const auto second = GetVertexPoint(mesh, vertices[1]);
    const auto [t0, d0] = SegmentCoordinateAndDistance(first, edge);
    const auto [t1, d1] = SegmentCoordinateAndDistance(second, edge);
    if (d0 > distance_tol || d1 > distance_tol)
    {
      continue;
    }
    const double t_min = std::min(t0, t1);
    const double t_max = std::max(t0, t1);
    if (t_min < -coordinate_tol || t_max > 1.0 + coordinate_tol)
    {
      continue;
    }
    const double sign =
        (Dot(Subtract(second, first), Subtract(edge[1], edge[0])) >= 0.0) ? 1.0 : -1.0;
    candidates.push_back(Candidate{edge_idx, std::max(0.0, t_min), std::min(1.0, t_max),
                                   Distance(first, second), sign});
  }
  std::sort(candidates.begin(), candidates.end(),
            [](const auto &a, const auto &b) { return a.t_min < b.t_min; });

  int local_count = static_cast<int>(candidates.size());
  int global_count = local_count;
  Mpi::GlobalSum(1, &global_count, mesh.GetComm());
  MFEM_VERIFY(global_count > 0,
              "\"TerminalEdges\" entry did not match any mesh edge chain!");

  double local_length = 0.0;
  for (const auto &candidate : candidates)
  {
    local_length += candidate.length;
  }
  double global_length = local_length;
  Mpi::GlobalSum(1, &global_length, mesh.GetComm());

  LumpedPortData::TerminalEdgeChain chain;
  chain.endpoints = edge;
  chain.edge_count = global_count;
  chain.length = global_length;
  chain.mesh_edges.reserve(candidates.size());
  chain.mesh_edge_signs.reserve(candidates.size());
  for (const auto &candidate : candidates)
  {
    chain.mesh_edges.push_back(candidate.edge);
    chain.mesh_edge_signs.push_back(candidate.sign);
  }
  return chain;
}

std::array<TerminalEdge, 2> BuildTerminalVoltageEdges(const TerminalEdge &first,
                                                      const TerminalEdge &second)
{
  const double direct = Distance(first[0], second[0]) + Distance(first[1], second[1]);
  const double crossed = Distance(first[0], second[1]) + Distance(first[1], second[0]);
  if (direct <= crossed)
  {
    return {TerminalEdge{first[0], second[0]}, TerminalEdge{first[1], second[1]}};
  }
  return {TerminalEdge{first[0], second[1]}, TerminalEdge{first[1], second[0]}};
}

void AddTerminalEdgeChainFunctional(const LumpedPortData::TerminalEdgeChain &chain,
                                    const mfem::ParFiniteElementSpace &nd_fespace,
                                    mfem::Vector &lf, double coeff)
{
  MFEM_VERIFY(chain.mesh_edges.size() == chain.mesh_edge_signs.size(),
              "Terminal edge chain metadata size mismatch!");
  mfem::Array<int> dofs;
  for (std::size_t i = 0; i < chain.mesh_edges.size(); i++)
  {
    nd_fespace.GetEdgeDofs(chain.mesh_edges[i], dofs);
    MFEM_VERIFY(dofs.Size() == 1,
                "\"TerminalEdges\" driven ports currently require first-order ND edge "
                "DOFs!");
    double dof_sign = 1.0;
    const int ldof = mfem::FiniteElementSpace::DecodeDof(dofs[0], dof_sign);
    lf(ldof) += coeff * chain.mesh_edge_signs[i] * dof_sign;
  }
}

double IntegrateTerminalEdgeChain(const LumpedPortData::TerminalEdgeChain &chain,
                                  const mfem::ParGridFunction &field)
{
  MFEM_VERIFY(chain.mesh_edges.size() == chain.mesh_edge_signs.size(),
              "Terminal edge chain metadata size mismatch!");
  const auto &nd_fespace = *field.ParFESpace();
  const double *values = field.HostRead();
  mfem::Array<int> dofs;
  double value = 0.0;
  for (std::size_t i = 0; i < chain.mesh_edges.size(); i++)
  {
    nd_fespace.GetEdgeDofs(chain.mesh_edges[i], dofs);
    MFEM_VERIFY(dofs.Size() == 1,
                "\"TerminalEdges\" driven ports currently require first-order ND edge "
                "DOFs!");
    double dof_sign = 1.0;
    const int ldof = mfem::FiniteElementSpace::DecodeDof(dofs[0], dof_sign);
    value += chain.mesh_edge_signs[i] * dof_sign * values[ldof];
  }
  return value;
}

bool VertexOnTerminalEdge(const Point &x, const TerminalEdge &edge)
{
  constexpr double distance_tol = 1.0e-8;
  constexpr double coordinate_tol = 1.0e-8;
  const auto [t, distance] = SegmentCoordinateAndDistance(x, edge);
  return distance <= distance_tol && t >= -coordinate_tol && t <= 1.0 + coordinate_tol;
}

class TerminalSheetModeCoefficient : public mfem::VectorCoefficient
{
private:
  const mfem::ParGridFunction &potential;
  const mfem::ParSubMesh &submesh;
  const std::unordered_map<int, int> &submesh_parent_elems;
  const std::unordered_set<int> &selected_parent_elems;
  mfem::IsoparametricTransformation T_loc;
  double scaling;

public:
  TerminalSheetModeCoefficient(const mfem::ParGridFunction &potential,
                               const mfem::ParSubMesh &submesh,
                               const std::unordered_map<int, int> &submesh_parent_elems,
                               const std::unordered_set<int> &selected_parent_elems,
                               double scaling = 1.0)
    : mfem::VectorCoefficient(submesh.SpaceDimension()), potential(potential),
      submesh(submesh), submesh_parent_elems(submesh_parent_elems),
      selected_parent_elems(selected_parent_elems), scaling(scaling)
  {
  }

  void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    mfem::ElementTransformation *T_submesh = nullptr;
    if (T.mesh == submesh.GetParent())
    {
      MFEM_ASSERT(T.ElementType == mfem::ElementTransformation::BDR_ELEMENT,
                  "TerminalSheetModeCoefficient requires ElementType::BDR_ELEMENT when "
                  "not used on a SubMesh!");
      auto it = submesh_parent_elems.find(T.ElementNo);
      if (it == submesh_parent_elems.end() ||
          selected_parent_elems.find(T.ElementNo) == selected_parent_elems.end())
      {
        V.SetSize(vdim);
        V = 0.0;
        return;
      }
      submesh.GetElementTransformation(it->second, &T_loc);
      T_loc.SetIntPoint(&ip);
      T_submesh = &T_loc;
    }
    else if (T.mesh == &submesh)
    {
      MFEM_ASSERT(T.ElementType == mfem::ElementTransformation::ELEMENT,
                  "TerminalSheetModeCoefficient requires ElementType::ELEMENT when used "
                  "on a SubMesh!");
      T_submesh = &T;
    }
    else
    {
      MFEM_ABORT("Invalid mesh for TerminalSheetModeCoefficient!");
    }

    potential.GetGradient(*T_submesh, V);
    V *= -scaling;
  }
};

class TerminalModalDampingOperator : public Operator
{
private:
  std::vector<Vector> modes;

public:
  TerminalModalDampingOperator(int size) : Operator(size) {}

  void AddMode(const Vector &mode)
  {
    modes.emplace_back(mode.Size());
    modes.back() = mode;
    modes.back().UseDevice(true);
  }

  void Mult(const Vector &x, Vector &y) const override
  {
    y = 0.0;
    AddMult(x, y);
  }

  void MultTranspose(const Vector &x, Vector &y) const override { Mult(x, y); }

  void AddMult(const Vector &x, Vector &y, const double a = 1.0) const override
  {
    for (const auto &mode : modes)
    {
      y.Add(a * (mode * x), mode);
    }
  }

  void AddMultTranspose(const Vector &x, Vector &y, const double a = 1.0) const override
  {
    AddMult(x, y, a);
  }

  void AssembleDiagonal(Vector &diag) const override
  {
    diag.SetSize(width);
    diag = 0.0;
    for (const auto &mode : modes)
    {
      const double *m = mode.HostRead();
      double *d = diag.HostReadWrite();
      for (int i = 0; i < mode.Size(); i++)
      {
        d[i] += m[i] * m[i];
      }
    }
  }
};

}  // namespace

struct LumpedPortData::TerminalSheetMode
{
  mfem::Array<int> attr_list;
  std::unique_ptr<Mesh> port_mesh;
  std::unique_ptr<mfem::FiniteElementCollection> port_h1_fec;
  std::unique_ptr<FiniteElementSpace> port_h1_fespace;
  std::unique_ptr<mfem::ParGridFunction> potential;
  std::unordered_map<int, int> submesh_parent_elems;
  std::unordered_set<int> selected_parent_elems;
  std::unordered_set<int> selected_submesh_elems;
  double norm_sq = 0.0;
  double selected_relative_permittivity = 0.0;
  int selected_adjacent_attr = 0;
  int global_selected_elements = 0;
  int global_signal_vertices = 0;
  int global_reference_vertices = 0;
  int global_ess_tdofs = 0;

  TerminalSheetMode(const mfem::Array<int> &attrs,
                    const std::array<TerminalEdge, 2> &terminals, const mfem::ParMesh &mesh,
                    const MaterialOperator &mat_op)
  {
    attr_list.Append(attrs);
    port_mesh = std::make_unique<Mesh>(std::make_unique<mfem::ParSubMesh>(
        mfem::ParSubMesh::CreateFromBoundary(mesh, attr_list)));
    port_h1_fec = std::make_unique<mfem::H1_FECollection>(1, port_mesh->Dimension());
    port_h1_fespace = std::make_unique<FiniteElementSpace>(*port_mesh, port_h1_fec.get());
    potential = std::make_unique<mfem::ParGridFunction>(&port_h1_fespace->Get());

    const auto &port_submesh = static_cast<const mfem::ParSubMesh &>(port_mesh->Get());
    const mfem::Array<int> &parent_elems = port_submesh.GetParentElementIDMap();
    for (int i = 0; i < parent_elems.Size(); i++)
    {
      submesh_parent_elems[parent_elems[i]] = i;
    }

    SelectAdjacentSide(mesh, mat_op);
    SolvePotential(terminals);
    norm_sq = ComputeNormSq();
    const auto [phi_min, phi_max] = ComputePotentialRange();
    Mpi::Print("\nTerminal sheet mode diagnostics:"
               " elements = {:d}, selected elements = {:d}, selected adjacent attr = {:d}, "
               "selected ε_r = {:.6e}, vertices = {:d}, signal vertices = {:d}, "
               "reference vertices = {:d}, essential true DOFs = {:d}, phi = "
               "[{:.6e}, {:.6e}], "
               "∫|E_1V|² dS = {:.6e}\n",
               port_mesh->GetNE(), global_selected_elements, selected_adjacent_attr,
               selected_relative_permittivity, port_mesh->Get().GetNV(),
               global_signal_vertices, global_reference_vertices, global_ess_tdofs, phi_min,
               phi_max, norm_sq);
    MFEM_VERIFY(norm_sq > 0.0, "Terminal sheet mode produced zero electric-field norm!");
  }

  static double RelativePermittivityScalar(const MaterialOperator &mat_op, int attr)
  {
    const auto eps = mat_op.GetPermittivityReal(attr);
    const int n = std::min(eps.Height(), eps.Width());
    double trace = 0.0;
    for (int i = 0; i < n; i++)
    {
      trace += eps(i, i);
    }
    return trace / static_cast<double>(n);
  }

  void SelectAdjacentSide(const mfem::ParMesh &parent_mesh, const MaterialOperator &mat_op)
  {
    int domain_attr_max = parent_mesh.attributes.Size() ? parent_mesh.attributes.Max() : 0;
    Mpi::GlobalMax(1, &domain_attr_max, parent_mesh.GetComm());

    std::vector<long long> adjacent_counts(domain_attr_max + 1, 0);
    std::vector<double> adjacent_eps(domain_attr_max + 1, mfem::infinity());
    std::unordered_map<int, int> parent_adjacent_attr;
    parent_adjacent_attr.reserve(submesh_parent_elems.size());

    for (const auto &[parent_be, submesh_elem] : submesh_parent_elems)
    {
      int elem_id = -1, info = 0;
      parent_mesh.GetBdrElementAdjacentElement(parent_be, elem_id, info);
      if (elem_id < 0)
      {
        continue;
      }
      const int domain_attr = parent_mesh.GetAttribute(elem_id);
      if (domain_attr <= 0 || domain_attr > domain_attr_max)
      {
        continue;
      }
      parent_adjacent_attr[parent_be] = domain_attr;
      adjacent_counts[domain_attr]++;
      adjacent_eps[domain_attr] = std::min(adjacent_eps[domain_attr],
                                           RelativePermittivityScalar(mat_op, domain_attr));
    }

    Mpi::GlobalSum(static_cast<int>(adjacent_counts.size()), adjacent_counts.data(),
                   parent_mesh.GetComm());
    Mpi::GlobalMin(static_cast<int>(adjacent_eps.size()), adjacent_eps.data(),
                   parent_mesh.GetComm());

    selected_adjacent_attr = 0;
    selected_relative_permittivity = mfem::infinity();
    for (int attr = 1; attr <= domain_attr_max; attr++)
    {
      if (adjacent_counts[attr] == 0)
      {
        continue;
      }
      if (adjacent_eps[attr] < selected_relative_permittivity ||
          (adjacent_eps[attr] == selected_relative_permittivity &&
           (selected_adjacent_attr == 0 || attr < selected_adjacent_attr)))
      {
        selected_adjacent_attr = attr;
        selected_relative_permittivity = adjacent_eps[attr];
      }
    }
    MFEM_VERIFY(selected_adjacent_attr > 0,
                "Terminal sheet mode did not find an adjacent material side!");

    for (const auto &[parent_be, domain_attr] : parent_adjacent_attr)
    {
      if (domain_attr != selected_adjacent_attr)
      {
        continue;
      }
      selected_parent_elems.insert(parent_be);
      selected_submesh_elems.insert(submesh_parent_elems.at(parent_be));
    }
    global_selected_elements = static_cast<int>(selected_parent_elems.size());
    Mpi::GlobalSum(1, &global_selected_elements, parent_mesh.GetComm());
    MFEM_VERIFY(global_selected_elements > 0,
                "Terminal sheet mode adjacent-side selection produced no elements!");
  }

  void SolvePotential(const std::array<TerminalEdge, 2> &terminals)
  {
    auto &fespace = port_h1_fespace->Get();
    const auto &mesh = *fespace.GetParMesh();
    MFEM_VERIFY(fespace.GetMaxElementOrder() == 1,
                "Terminal sheet mode currently requires first-order H1 elements!");
    *potential = 0.0;

    mfem::Array<int> dofs;
    mfem::Array<int> vertices;
    mfem::Array<int> edges;
    mfem::Array<int> orientations;
    std::unordered_map<int, int> selected_edge_counts;
    for (int elem : selected_submesh_elems)
    {
      mesh.GetElementEdges(elem, edges, orientations);
      for (int i = 0; i < edges.Size(); i++)
      {
        selected_edge_counts[edges[i]]++;
      }
    }
    std::unordered_set<int> boundary_vertices;
    for (const auto &[edge, count] : selected_edge_counts)
    {
      if (count != 1)
      {
        continue;
      }
      mesh.GetEdgeVertices(edge, vertices);
      for (int i = 0; i < vertices.Size(); i++)
      {
        boundary_vertices.insert(vertices[i]);
      }
    }

    std::set<int> ess_tdofs;
    std::set<int> signal_tdofs;
    int local_signal_vertices = 0, local_reference_vertices = 0;
    for (int v = 0; v < mesh.GetNV(); v++)
    {
      const auto point = GetVertexPoint(mesh, v);
      const bool on_signal = VertexOnTerminalEdge(point, terminals[0]);
      const bool on_reference =
          !on_signal && boundary_vertices.find(v) != boundary_vertices.end();
      MFEM_VERIFY(!(on_signal && on_reference),
                  "\"TerminalEdges\" signal and reference chains overlap on the port "
                  "sheet!");
      if (!on_signal && !on_reference)
      {
        continue;
      }
      local_signal_vertices += on_signal ? 1 : 0;
      local_reference_vertices += on_reference ? 1 : 0;
      const double value = on_signal ? 1.0 : 0.0;
      fespace.GetVertexDofs(v, dofs);
      for (int i = 0; i < dofs.Size(); i++)
      {
        double sign = 1.0;
        const int ldof = mfem::FiniteElementSpace::DecodeDof(dofs[i], sign);
        (*potential)(ldof) = sign * value;
        const int ltdof = fespace.GetLocalTDofNumber(ldof);
        if (ltdof >= 0)
        {
          ess_tdofs.insert(ltdof);
          if (on_signal)
          {
            signal_tdofs.insert(ltdof);
          }
        }
      }
    }

    int global_counts[2] = {local_signal_vertices, local_reference_vertices};
    Mpi::GlobalSum(2, global_counts, mesh.GetComm());
    global_signal_vertices = global_counts[0];
    global_reference_vertices = global_counts[1];
    MFEM_VERIFY(global_counts[0] > 0 && global_counts[1] > 0,
                "\"TerminalEdges\" did not map to both signal and reference vertices on "
                "the port sheet!");

    mfem::Array<int> ess_tdof_list;
    ess_tdof_list.Reserve(static_cast<int>(ess_tdofs.size()));
    for (int tdof : ess_tdofs)
    {
      ess_tdof_list.Append(tdof);
    }
    int local_ess_tdofs = ess_tdof_list.Size();
    global_ess_tdofs = local_ess_tdofs;
    Mpi::GlobalSum(1, &global_ess_tdofs, mesh.GetComm());
    MFEM_VERIFY(global_ess_tdofs > 0,
                "\"TerminalEdges\" did not map to owned H1 true DOFs on the port sheet!");

    mfem::Vector true_potential(fespace.GetTrueVSize());
    true_potential.UseDevice(true);
    true_potential = 0.0;
    {
      double *values = true_potential.HostReadWrite();
      for (int tdof : signal_tdofs)
      {
        values[tdof] = 1.0;
      }
    }
    potential->SetFromTrueDofs(true_potential);

    mfem::ConstantCoefficient one(1.0);
    mfem::ParBilinearForm a(&fespace);
    a.AddDomainIntegrator(new mfem::DiffusionIntegrator(one));
    a.Assemble();
    a.Finalize();

    mfem::ParLinearForm b(&fespace);
    b = 0.0;

    mfem::HypreParMatrix A;
    mfem::Vector X, B;
    a.FormLinearSystem(ess_tdof_list, *potential, b, A, X, B);

    mfem::CGSolver cg(mesh.GetComm());
    cg.SetPrintLevel(0);
    cg.SetRelTol(1.0e-12);
    cg.SetAbsTol(1.0e-14);
    cg.SetMaxIter(500);
    mfem::HypreBoomerAMG amg(A);
    amg.SetPrintLevel(0);
    cg.SetPreconditioner(amg);
    cg.SetOperator(A);
    cg.Mult(B, X);
    a.RecoverFEMSolution(X, b, *potential);
  }

  double ComputeNormSq()
  {
    auto &mesh = port_mesh->Get();
    double local_norm = 0.0;
    mfem::Vector grad(mesh.SpaceDimension());
    const int order = 2 * port_h1_fespace->GetMaxElementOrder() + 2;
    for (int i = 0; i < mesh.GetNE(); i++)
    {
      if (selected_submesh_elems.find(i) == selected_submesh_elems.end())
      {
        continue;
      }
      auto *T = mesh.GetElementTransformation(i);
      const mfem::IntegrationRule &ir = mfem::IntRules.Get(T->GetGeometryType(), order);
      for (int j = 0; j < ir.GetNPoints(); j++)
      {
        const mfem::IntegrationPoint &ip = ir.IntPoint(j);
        T->SetIntPoint(&ip);
        potential->GetGradient(*T, grad);
        local_norm += ip.weight * T->Weight() * (grad * grad);
      }
    }
    double global_norm = local_norm;
    Mpi::GlobalSum(1, &global_norm, mesh.GetComm());
    return global_norm;
  }

  std::pair<double, double> ComputePotentialRange() const
  {
    double phi_min = mfem::infinity();
    double phi_max = -mfem::infinity();
    const double *values = potential->HostRead();
    for (int i = 0; i < potential->Size(); i++)
    {
      phi_min = std::min(phi_min, values[i]);
      phi_max = std::max(phi_max, values[i]);
    }
    Mpi::GlobalMin(1, &phi_min, potential->ParFESpace()->GetComm());
    Mpi::GlobalMax(1, &phi_max, potential->ParFESpace()->GetComm());
    return {phi_min, phi_max};
  }

  std::unique_ptr<mfem::VectorCoefficient> GetCoefficient(double coeff) const
  {
    const auto &port_submesh = static_cast<const mfem::ParSubMesh &>(port_mesh->Get());
    return std::make_unique<RestrictedVectorCoefficient<TerminalSheetModeCoefficient>>(
        attr_list, *potential, port_submesh, submesh_parent_elems, selected_parent_elems,
        coeff);
  }
};

LumpedPortData::~LumpedPortData() = default;

LumpedPortData::LumpedPortData(const config::LumpedPortData &data,
                               const MaterialOperator &mat_op, const mfem::ParMesh &mesh)
  : mat_op(mat_op), excitation(data.excitation), active(data.active)
{
  // Check inputs. Only one of the circuit or per square properties should be specified
  // for the port boundary.
  bool has_circ = (std::abs(data.R) + std::abs(data.L) + std::abs(data.C) > 0.0);
  bool has_surf = (std::abs(data.Rs) + std::abs(data.Ls) + std::abs(data.Cs) > 0.0);
  MFEM_VERIFY(has_circ || has_surf,
              "Lumped port boundary has no R/L/C or Rs/Ls/Cs defined, needs "
              "at least one!");
  MFEM_VERIFY(!(has_circ && has_surf),
              "Lumped port boundary has both R/L/C and Rs/Ls/Cs defined, "
              "should only use one!");

  if (HasExcitation())
  {
    if (has_circ)
    {
      MFEM_VERIFY(data.R > 0.0, "Excited lumped port must have nonzero resistance!");
      MFEM_VERIFY(data.C == 0.0 && data.L == 0.0,
                  "Lumped port excitations do not support nonzero reactance!");
    }
    else
    {
      MFEM_VERIFY(data.Rs > 0.0, "Excited lumped port must have nonzero resistance!");
      MFEM_VERIFY(data.Cs == 0.0 && data.Ls == 0.0,
                  "Lumped port excitations do not support nonzero reactance!");
    }
  }

  // Construct the port elements allowing for a possible multielement lumped port.
  for (const auto &elem : data.elements)
  {
    mfem::Array<int> attr_list;
    attr_list.Append(elem.attributes.data(), elem.attributes.size());
    switch (elem.coordinate_system)
    {
      case CoordinateSystem::CYLINDRICAL:
        elems.push_back(
            std::make_unique<CoaxialElementData>(elem.direction, attr_list, mesh));
        break;
      case CoordinateSystem::CARTESIAN:
        elems.push_back(std::make_unique<UniformElementData>(
            elem.direction, attr_list, mesh, elem.length, elem.width));
        break;
    }
    if (!elem.terminal_edges.empty())
    {
      MFEM_VERIFY(elem.terminal_edges.size() == 2,
                  "\"TerminalEdges\" must contain exactly two endpoint pairs!");
      PrintBoundaryAdjacencyDiagnostics(
          mesh, attr_list, fmt::format("terminal lumped port element {:d}", elems.size()));
      terminal_edges.push_back({FindTerminalEdgeChain(mesh, elem.terminal_edges[0]),
                                FindTerminalEdgeChain(mesh, elem.terminal_edges[1])});
      const auto voltage_edges =
          BuildTerminalVoltageEdges(elem.terminal_edges[0], elem.terminal_edges[1]);
      terminal_voltage_edges.push_back({FindTerminalEdgeChain(mesh, voltage_edges[0]),
                                        FindTerminalEdgeChain(mesh, voltage_edges[1])});
      terminal_sheet_modes.push_back(std::make_unique<TerminalSheetMode>(
          attr_list,
          std::array<TerminalEdge, 2>{elem.terminal_edges[0], elem.terminal_edges[1]}, mesh,
          mat_op));
    }
    else
    {
      terminal_sheet_modes.push_back(nullptr);
    }
  }

  // Populate the property data for the lumped port.
  if (std::abs(data.Rs) + std::abs(data.Ls) + std::abs(data.Cs) == 0.0)
  {
    R = data.R;
    L = data.L;
    C = data.C;
  }
  else
  {
    // If defined by surface properties, need to compute circuit properties for the
    // multielement port.
    double ooR = 0.0, ooL = 0.0;
    R = L = C = 0.0;
    for (const auto &elem : elems)
    {
      const double sq = elem->GetGeometryWidth() / elem->GetGeometryLength();
      if (std::abs(data.Rs) > 0.0)
      {
        ooR += sq / data.Rs;
      }
      if (std::abs(data.Ls) > 0.0)
      {
        ooL += sq / data.Ls;
      }
      if (std::abs(data.Cs) > 0.0)
      {
        C += sq * data.Cs;
      }
    }
    if (std::abs(ooR) > 0.0)
    {
      R = 1.0 / ooR;
    }
    if (std::abs(ooL) > 0.0)
    {
      L = 1.0 / ooL;
    }
  }

  if (HasTerminalEdges())
  {
    fmt::memory_buffer buffer{};
    auto out = fmt::appender{buffer};
    fmt::format_to(out, "\nResolved HFSS-style terminal edge chains for lumped port:\n");
    for (std::size_t element_idx = 0; element_idx < terminal_edges.size(); element_idx++)
    {
      const auto &chains = terminal_edges[element_idx];
      fmt::format_to(out,
                     " Element {:d}: signal edges = {:d} (length {:.6e}), reference "
                     "edges = {:d} (length {:.6e})\n",
                     element_idx + 1, chains[0].edge_count, chains[0].length,
                     chains[1].edge_count, chains[1].length);
    }
    fmt::format_to(out, "Resolved terminal voltage edge chains for lumped port:\n");
    for (std::size_t element_idx = 0; element_idx < terminal_voltage_edges.size();
         element_idx++)
    {
      const auto &chains = terminal_voltage_edges[element_idx];
      fmt::format_to(out,
                     " Element {:d}: voltage path A edges = {:d} (length {:.6e}), "
                     "voltage path B edges = {:d} (length {:.6e})\n",
                     element_idx + 1, chains[0].edge_count, chains[0].length,
                     chains[1].edge_count, chains[1].length);
    }
    fmt::format_to(out, "Resolved terminal sheet modes for lumped port:\n");
    std::size_t terminal_mode_idx = 0;
    for (std::size_t element_idx = 0; element_idx < terminal_sheet_modes.size();
         element_idx++)
    {
      if (!terminal_sheet_modes[element_idx])
      {
        continue;
      }
      fmt::format_to(out,
                     " Element {:d}: ∫|E_1V|² dS = {:.6e}, effective squares = "
                     "{:.6e}\n",
                     ++terminal_mode_idx, terminal_sheet_modes[element_idx]->norm_sq,
                     GetToSquare(*elems[element_idx]));
    }
    Mpi::Print("{}", fmt::to_string(buffer));
  }
}

double LumpedPortData::GetToSquare(const LumpedElementData &elem) const
{
  for (std::size_t i = 0; i < elems.size(); i++)
  {
    if (elems[i].get() == &elem && i < terminal_sheet_modes.size() &&
        terminal_sheet_modes[i])
    {
      return terminal_sheet_modes[i]->norm_sq * elems.size();
    }
  }
  return elem.GetGeometryWidth() / elem.GetGeometryLength() * elems.size();
}

std::unique_ptr<mfem::VectorCoefficient>
LumpedPortData::GetModeCoefficient(std::size_t elem_idx, double coeff) const
{
  MFEM_VERIFY(elem_idx < elems.size(), "Invalid lumped port element index!");
  if (elem_idx < terminal_sheet_modes.size() && terminal_sheet_modes[elem_idx])
  {
    return terminal_sheet_modes[elem_idx]->GetCoefficient(coeff);
  }
  return elems[elem_idx]->GetModeCoefficient(coeff);
}

void AssemblePortModeLinearForm(const LumpedPortData &data,
                                mfem::ParFiniteElementSpace &nd_fespace, Vector &mode)
{
  const auto &mesh = *nd_fespace.GetParMesh();
  SumVectorCoefficient fb(mesh.SpaceDimension());
  mfem::Array<int> attr_list;
  for (std::size_t elem_idx = 0; elem_idx < data.elems.size(); elem_idx++)
  {
    const auto &elem = *data.elems[elem_idx];
    const double Rs = data.R * data.GetToSquare(elem);
    const bool terminal_mode =
        elem_idx < data.terminal_sheet_modes.size() && data.terminal_sheet_modes[elem_idx];
    const double Hinc =
        (std::abs(Rs) > 0.0)
            ? (terminal_mode
                   ? std::sqrt(data.R) / Rs
                   : 1.0 / std::sqrt(Rs * elem.GetGeometryWidth() *
                                     elem.GetGeometryLength() * data.elems.size()))
            : 0.0;
    fb.AddCoefficient(data.GetModeCoefficient(elem_idx, Hinc));
    attr_list.Append(elem.GetAttrList());
  }

  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> attr_marker = mesh::AttrToMarker(bdr_attr_max, attr_list);
  mfem::LinearForm lf(&nd_fespace);
  lf.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fb), attr_marker);
  lf.UseFastAssembly(false);
  lf.UseDevice(false);
  lf.Assemble();
  lf.UseDevice(true);

  mode.SetSize(nd_fespace.GetVSize());
  mode.UseDevice(true);
  mode = lf;
}

std::complex<double>
LumpedPortData::GetCharacteristicImpedance(double omega,
                                           LumpedPortData::Branch branch) const
{
  MFEM_VERIFY((L == 0.0 && C == 0.0) || branch == Branch::R || omega > 0.0,
              "Lumped port with nonzero reactance requires frequency in order to define "
              "characteristic impedance!");
  std::complex<double> Y = 0.0;
  if (std::abs(R) > 0.0 && (branch == Branch::TOTAL || branch == Branch::R))
  {
    Y += 1.0 / R;
  }
  if (std::abs(L) > 0.0 && (branch == Branch::TOTAL || branch == Branch::L))
  {
    Y += 1.0 / (1i * omega * L);
  }
  if (std::abs(C) > 0.0 && (branch == Branch::TOTAL || branch == Branch::C))
  {
    Y += 1i * omega * C;
  }
  MFEM_VERIFY(std::abs(Y) > 0.0,
              "Characteristic impedance requested for lumped port with zero admittance!")
  return 1.0 / Y;
}

double LumpedPortData::GetExcitationPower() const
{
  // The lumped port excitation is normalized such that the power integrated over the port
  // is 1: ∫ (E_inc x H_inc) ⋅ n dS = 1.
  return HasExcitation() ? 1.0 : 0.0;
}

double LumpedPortData::GetExcitationVoltage() const
{
  // Incident voltage should be the same across all elements of an excited lumped port.
  if (HasExcitation())
  {
    if (HasTerminalEdges())
    {
      return std::sqrt(R);
    }
    double V_inc = 0.0;
    for (const auto &elem : elems)
    {
      const double Rs = R * GetToSquare(*elem);
      const double E_inc = std::sqrt(
          Rs / (elem->GetGeometryWidth() * elem->GetGeometryLength() * elems.size()));
      V_inc += E_inc * elem->GetGeometryLength() / elems.size();
    }
    return V_inc;
  }
  else
  {
    return 0.0;
  }
}

void LumpedPortData::AddTerminalEdgeVoltageFunctional(
    const mfem::ParFiniteElementSpace &nd_fespace, Vector &lf, double coeff) const
{
  MFEM_VERIFY(HasTerminalEdges(),
              "Terminal edge voltage functional requested for a non-terminal port!");
  const double weight = 1.0 / (2.0 * static_cast<double>(terminal_voltage_edges.size()));
  for (const auto &edge_pair : terminal_voltage_edges)
  {
    AddTerminalEdgeChainFunctional(edge_pair[0], nd_fespace, lf, coeff * weight);
    AddTerminalEdgeChainFunctional(edge_pair[1], nd_fespace, lf, coeff * weight);
  }
}

void LumpedPortData::AddTerminalEdgeExcitationFunctional(
    const mfem::ParFiniteElementSpace &nd_fespace, Vector &lf, double coeff) const
{
  MFEM_VERIFY(HasTerminalEdges(),
              "Terminal edge excitation requested for a non-terminal port!");
  const double weight = 1.0 / (2.0 * static_cast<double>(terminal_edges.size()));
  for (const auto &edge_pair : terminal_edges)
  {
    AddTerminalEdgeChainFunctional(edge_pair[0], nd_fespace, lf, coeff * weight);
    AddTerminalEdgeChainFunctional(edge_pair[1], nd_fespace, lf, -coeff * weight);
  }
}

void LumpedPortData::InitializeLinearForms(mfem::ParFiniteElementSpace &nd_fespace) const
{
  const auto &mesh = *nd_fespace.GetParMesh();
  mfem::Array<int> attr_marker;
  if (!s || !v)
  {
    mfem::Array<int> attr_list;
    for (const auto &elem : elems)
    {
      attr_list.Append(elem->GetAttrList());
    }
    int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
    mesh::AttrToMarker(bdr_attr_max, attr_list, attr_marker);
  }

  // The port S-parameter, or the projection of the field onto the port mode, is computed
  // as: (E x H_inc) ⋅ n = E ⋅ (E_inc / Z_s), integrated over the port surface.
  if (!s)
  {
    SumVectorCoefficient fb(mesh.SpaceDimension());
    for (std::size_t elem_idx = 0; elem_idx < elems.size(); elem_idx++)
    {
      const auto &elem = *elems[elem_idx];
      const double Rs = R * GetToSquare(elem);
      const bool terminal_mode =
          elem_idx < terminal_sheet_modes.size() && terminal_sheet_modes[elem_idx];
      const double Hinc =
          (std::abs(Rs) > 0.0)
              ? (terminal_mode ? std::sqrt(R) / Rs
                               : 1.0 / std::sqrt(Rs * elem.GetGeometryWidth() *
                                                 elem.GetGeometryLength() * elems.size()))
              : 0.0;
      fb.AddCoefficient(GetModeCoefficient(elem_idx, Hinc));
    }
    s = std::make_unique<mfem::LinearForm>(&nd_fespace);
    s->AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fb), attr_marker);
    s->UseFastAssembly(false);
    s->UseDevice(false);
    s->Assemble();
    s->UseDevice(true);
  }

  // The voltage across a port is computed using the electric field solution.
  // We have:
  //             V = ∫ E ⋅ l̂ dl = 1/w ∫ E ⋅ l̂ dS  (for rectangular ports)
  // or,
  //             V = 1/(2π) ∫ E ⋅ r̂ / r dS        (for coaxial ports).
  // We compute the surface integral via an inner product between the linear form with the
  // averaging function as a vector coefficient and the solution expansion coefficients.
  if (!v)
  {
    SumVectorCoefficient fb(mesh.SpaceDimension());
    for (const auto &elem : elems)
    {
      fb.AddCoefficient(
          elem->GetModeCoefficient(1.0 / (elem->GetGeometryWidth() * elems.size())));
    }
    v = std::make_unique<mfem::LinearForm>(&nd_fespace);
    v->AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fb), attr_marker);
    v->UseFastAssembly(false);
    v->UseDevice(false);
    v->Assemble();
    v->UseDevice(true);
  }
}

std::complex<double> LumpedPortData::GetPower(GridFunction &E, GridFunction &B) const
{
  // Compute port power, (E x H) ⋅ n = E ⋅ (-n x H), integrated over the port surface using
  // the computed E and H = μ⁻¹ B fields, where +n is the direction of propagation (into the
  // domain). The BdrSurfaceCurrentVectorCoefficient computes -n x H for an outward normal,
  // so we multiply by -1. The linear form is reconstructed from scratch each time due to
  // changing H.
  MFEM_VERIFY((E.HasImag() && B.HasImag()) || (!E.HasImag() && !B.HasImag()),
              "Mismatch between real- and complex-valued E and B fields in port power "
              "calculation!");
  const bool has_imag = E.HasImag();
  auto &nd_fespace = *E.ParFESpace();
  const auto &mesh = *nd_fespace.GetParMesh();
  SumVectorCoefficient fbr(mesh.SpaceDimension()), fbi(mesh.SpaceDimension());
  mfem::Array<int> attr_list;
  for (const auto &elem : elems)
  {
    fbr.AddCoefficient(
        std::make_unique<RestrictedVectorCoefficient<BdrSurfaceCurrentVectorCoefficient>>(
            elem->GetAttrList(), B.Real(), mat_op));
    if (has_imag)
    {
      fbi.AddCoefficient(
          std::make_unique<RestrictedVectorCoefficient<BdrSurfaceCurrentVectorCoefficient>>(
              elem->GetAttrList(), B.Imag(), mat_op));
    }
    attr_list.Append(elem->GetAttrList());
  }
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> attr_marker = mesh::AttrToMarker(bdr_attr_max, attr_list);
  std::complex<double> dot;
  {
    mfem::LinearForm pr(&nd_fespace);
    pr.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fbr), attr_marker);
    pr.UseFastAssembly(false);
    pr.UseDevice(false);
    pr.Assemble();
    pr.UseDevice(true);
    dot = -(pr * E.Real()) + (has_imag ? -1i * (pr * E.Imag()) : 0.0);
  }
  if (has_imag)
  {
    mfem::LinearForm pi(&nd_fespace);
    pi.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fbi), attr_marker);
    pi.UseFastAssembly(false);
    pi.UseDevice(false);
    pi.Assemble();
    pi.UseDevice(true);
    dot += -(pi * E.Imag()) + 1i * (pi * E.Real());
    Mpi::GlobalSum(1, &dot, E.ParFESpace()->GetComm());
    return dot;
  }
  else
  {
    double rdot = dot.real();
    Mpi::GlobalSum(1, &rdot, E.ParFESpace()->GetComm());
    return rdot;
  }
}

std::complex<double> LumpedPortData::GetSParameter(GridFunction &E) const
{
  // Compute port S-parameter, or the projection of the field onto the port mode.
  if (HasTerminalEdges())
  {
    InitializeLinearForms(*E.ParFESpace());
    std::complex<double> dot((*s) * E.Real(), 0.0);
    if (E.HasImag())
    {
      dot.imag((*s) * E.Imag());
    }
    Mpi::GlobalSum(1, &dot, E.GetComm());
    return dot;
  }
  InitializeLinearForms(*E.ParFESpace());
  std::complex<double> dot((*s) * E.Real(), 0.0);
  if (E.HasImag())
  {
    dot.imag((*s) * E.Imag());
  }
  Mpi::GlobalSum(1, &dot, E.GetComm());
  return dot;
}

std::complex<double> LumpedPortData::GetVoltage(GridFunction &E) const
{
  // Compute the average voltage across the port.
  if (HasTerminalEdges())
  {
    const double weight = 1.0 / (2.0 * static_cast<double>(terminal_voltage_edges.size()));
    std::complex<double> dot = 0.0;
    for (const auto &edge_pair : terminal_voltage_edges)
    {
      dot.real(dot.real() + weight * IntegrateTerminalEdgeChain(edge_pair[0], E.Real()));
      dot.real(dot.real() + weight * IntegrateTerminalEdgeChain(edge_pair[1], E.Real()));
      if (E.HasImag())
      {
        dot.imag(dot.imag() + weight * IntegrateTerminalEdgeChain(edge_pair[0], E.Imag()));
        dot.imag(dot.imag() + weight * IntegrateTerminalEdgeChain(edge_pair[1], E.Imag()));
      }
    }
    Mpi::GlobalSum(1, &dot, E.GetComm());
    return dot;
  }
  InitializeLinearForms(*E.ParFESpace());
  std::complex<double> dot((*v) * E.Real(), 0.0);
  if (E.HasImag())
  {
    dot.imag((*v) * E.Imag());
  }
  Mpi::GlobalSum(1, &dot, E.GetComm());
  return dot;
}

LumpedPortOperator::LumpedPortOperator(
    const std::map<int, config::LumpedPortData> &lumpedport, const Units &units,
    const MaterialOperator &mat_op, const mfem::ParMesh &mesh)
{
  SetUpBoundaryProperties(lumpedport, mat_op, mesh);
  PrintBoundaryInfo(units, mesh);
}

LumpedPortOperator::LumpedPortOperator(const IoData &iodata, const MaterialOperator &mat_op,
                                       const mfem::ParMesh &mesh)
  : LumpedPortOperator(iodata.boundaries.lumpedport, iodata.units, mat_op, mesh)
{
}

void LumpedPortOperator::SetUpBoundaryProperties(
    const std::map<int, config::LumpedPortData> &lumpedport, const MaterialOperator &mat_op,
    const mfem::ParMesh &mesh)
{
  // Check that lumped port boundary attributes have been specified correctly.
  if (!lumpedport.empty())
  {
    int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
    mfem::Array<int> bdr_attr_marker(bdr_attr_max), port_marker(bdr_attr_max);
    bdr_attr_marker = 0;
    port_marker = 0;
    for (auto attr : mesh.bdr_attributes)
    {
      bdr_attr_marker[attr - 1] = 1;
    }
    for (const auto &[idx, data] : lumpedport)
    {
      for (const auto &elem : data.elements)
      {
        for (auto attr : elem.attributes)
        {
          MFEM_VERIFY(attr > 0 && attr <= bdr_attr_max,
                      "Port boundary attribute tags must be non-negative and correspond to "
                      "boundaries in the mesh!");
          MFEM_VERIFY(bdr_attr_marker[attr - 1],
                      "Unknown port boundary attribute " << attr << "!");
          MFEM_VERIFY(!data.active || !port_marker[attr - 1],
                      "Boundary attribute is assigned to more than one lumped port!");
          port_marker[attr - 1] = 1;
        }
      }
    }
  }

  // Set up lumped port data structures.
  for (const auto &[idx, data] : lumpedport)
  {
    ports.try_emplace(idx, data, mat_op, mesh);
  }
}

void LumpedPortOperator::PrintBoundaryInfo(const Units &units, const mfem::ParMesh &mesh)
{
  if (ports.empty())
  {
    return;
  }
  fmt::memory_buffer buffer{};
  auto out = fmt::appender{buffer};
  using VT = Units::ValueType;

  // Print out BC info for all port attributes, for both active and inactive ports.
  fmt::format_to(out, "\nConfiguring Robin impedance BC for lumped ports at attributes:\n");
  for (const auto &[idx, data] : ports)
  {
    for (const auto &elem : data.elems)
    {
      for (auto attr : elem->GetAttrList())
      {
        fmt::format_to(out, " {:d}:", attr);
        if (std::abs(data.R) > 0.0)
        {
          double Rs = data.R * data.GetToSquare(*elem);
          fmt::format_to(out, " Rs = {:.3e} Ω/sq,",
                         units.Dimensionalize<VT::IMPEDANCE>(Rs));
        }
        if (std::abs(data.L) > 0.0)
        {
          double Ls = data.L * data.GetToSquare(*elem);
          fmt::format_to(out, " Ls = {:.3e} H/sq,",
                         units.Dimensionalize<VT::INDUCTANCE>(Ls));
        }
        if (std::abs(data.C) > 0.0)
        {
          double Cs = data.C / data.GetToSquare(*elem);
          fmt::format_to(out, " Cs = {:.3e} F/sq,",
                         units.Dimensionalize<VT::CAPACITANCE>(Cs));
        }
        fmt::format_to(out, " n = ({:+.1f})\n",
                       fmt::join(mesh::GetSurfaceNormal(mesh, attr), ","));
      }
    }
  }

  // Print out port info for all active ports.
  fmt::memory_buffer buffer_active{};
  for (const auto &[idx, data] : ports)
  {
    if (!data.active)
    {
      continue;
    }
    std::vector<std::string> active_port_str;
    if (std::abs(data.R) > 0.0)
    {
      active_port_str.emplace_back(
          fmt::format("R = {:.3e} Ω", units.Dimensionalize<VT::IMPEDANCE>(data.R)));
    }
    if (std::abs(data.L) > 0.0)
    {
      active_port_str.emplace_back(
          fmt::format("L = {:.3e} H", units.Dimensionalize<VT::INDUCTANCE>(data.L)));
    }
    if (std::abs(data.C) > 0.0)
    {
      active_port_str.emplace_back(
          fmt::format("C = {:.3e} F", units.Dimensionalize<VT::CAPACITANCE>(data.C)));
    }
    fmt::format_to(fmt::appender{buffer_active}, " Index = {:d}: {}\n", idx,
                   fmt::join(active_port_str, ", "));
  }
  if (buffer_active.size() > 0)
  {
    fmt::format_to(out, "\nConfiguring lumped port circuit properties:\n");
    buffer.append(buffer_active);
    buffer_active.clear();
  }

  // Print some information for excited lumped ports.
  for (const auto &[idx, data] : ports)
  {
    if (!data.HasExcitation())
    {
      continue;
    }

    for (const auto &elem : data.elems)
    {
      for (auto attr : elem->GetAttrList())
      {
        fmt::format_to(fmt::appender{buffer_active}, " {:d}: Index = {:d}\n", attr, idx);
      }
    }
  }
  if (buffer_active.size() > 0)
  {
    fmt::format_to(out,
                   "\nConfiguring lumped port excitation source term at attributes:\n");
    buffer.append(buffer_active);
  }

  Mpi::Print("{}", fmt::to_string(buffer));
}

const LumpedPortData &LumpedPortOperator::GetPort(int idx) const
{
  auto it = ports.find(idx);
  MFEM_VERIFY(it != ports.end(), "Unknown lumped port index requested!");
  return it->second;
}

mfem::Array<int> LumpedPortOperator::GetAttrList() const
{
  mfem::Array<int> attr_list;
  for (const auto &[idx, data] : ports)
  {
    if (!data.active)
    {
      continue;
    }
    for (const auto &elem : data.elems)
    {
      attr_list.Append(elem->GetAttrList());
    }
  }
  return attr_list;
}

mfem::Array<int> LumpedPortOperator::GetRsAttrList() const
{
  mfem::Array<int> attr_list;
  for (const auto &[idx, data] : ports)
  {
    if (!data.active)
    {
      continue;
    }
    if (data.HasTerminalEdges())
    {
      continue;
    }
    if (std::abs(data.R) > 0.0)
    {
      for (const auto &elem : data.elems)
      {
        attr_list.Append(elem->GetAttrList());
      }
    }
  }
  return attr_list;
}

mfem::Array<int> LumpedPortOperator::GetLsAttrList() const
{
  mfem::Array<int> attr_list;
  for (const auto &[idx, data] : ports)
  {
    if (!data.active)
    {
      continue;
    }
    if (std::abs(data.L) > 0.0)
    {
      for (const auto &elem : data.elems)
      {
        attr_list.Append(elem->GetAttrList());
      }
    }
  }
  return attr_list;
}

mfem::Array<int> LumpedPortOperator::GetCsAttrList() const
{
  mfem::Array<int> attr_list;
  for (const auto &[idx, data] : ports)
  {
    if (!data.active)
    {
      continue;
    }
    if (std::abs(data.C) > 0.0)
    {
      for (const auto &elem : data.elems)
      {
        attr_list.Append(elem->GetAttrList());
      }
    }
  }
  return attr_list;
}

void LumpedPortOperator::AddStiffnessBdrCoefficients(double coeff,
                                                     MaterialPropertyCoefficient &fb)
{
  // Add lumped inductor boundaries to the bilinear form.
  for (const auto &[idx, data] : ports)
  {
    if (!data.active)
    {
      continue;
    }
    if (std::abs(data.L) > 0.0)
    {
      for (const auto &elem : data.elems)
      {
        const double Ls = data.L * data.GetToSquare(*elem);
        fb.AddMaterialProperty(data.mat_op.GetCeedBdrAttributes(elem->GetAttrList()),
                               coeff / Ls);
      }
    }
  }
}

void LumpedPortOperator::AddDampingBdrCoefficients(double coeff,
                                                   MaterialPropertyCoefficient &fb)
{
  // Add lumped resistor boundaries to the bilinear form.
  for (const auto &[idx, data] : ports)
  {
    if (!data.active)
    {
      continue;
    }
    if (data.HasTerminalEdges())
    {
      continue;
    }
    if (std::abs(data.R) > 0.0)
    {
      for (const auto &elem : data.elems)
      {
        const double Rs = data.R * data.GetToSquare(*elem);
        fb.AddMaterialProperty(data.mat_op.GetCeedBdrAttributes(elem->GetAttrList()),
                               coeff / Rs);
      }
    }
  }
}

std::unique_ptr<Operator> LumpedPortOperator::GetTerminalModalDampingOperator(
    mfem::ParFiniteElementSpace &nd_fespace) const
{
  std::unique_ptr<TerminalModalDampingOperator> op;
  Vector mode;
  for (const auto &[idx, data] : ports)
  {
    if (!data.active || !data.HasTerminalEdges() || std::abs(data.R) == 0.0)
    {
      continue;
    }
    AssemblePortModeLinearForm(data, nd_fespace, mode);
    if (!op)
    {
      op = std::make_unique<TerminalModalDampingOperator>(nd_fespace.GetVSize());
    }
    op->AddMode(mode);
  }
  return op;
}

void LumpedPortOperator::AddMassBdrCoefficients(double coeff,
                                                MaterialPropertyCoefficient &fb)
{
  // Add lumped capacitance boundaries to the bilinear form.
  for (const auto &[idx, data] : ports)
  {
    if (!data.active)
    {
      continue;
    }
    if (std::abs(data.C) > 0.0)
    {
      for (const auto &elem : data.elems)
      {
        const double Cs = data.C / data.GetToSquare(*elem);
        fb.AddMaterialProperty(data.mat_op.GetCeedBdrAttributes(elem->GetAttrList()),
                               coeff * Cs);
      }
    }
  }
}

void LumpedPortOperator::AddExcitationBdrCoefficients(int excitation_idx,
                                                      SumVectorCoefficient &fb)
{
  // Construct the RHS source term for lumped port boundaries, which looks like -U_inc =
  // +2 iω/Z_s E_inc for a port boundary with an incident field E_inc. The chosen incident
  // field magnitude corresponds to a unit incident power over the full port boundary. See
  // p. 49 and p. 82 of the COMSOL RF Module manual for more detail.
  // Note: The real RHS returned here does not yet have the factor of (iω) included, so
  // works for time domain simulations requiring RHS -U_inc(t).
  for (const auto &[idx, data] : ports)
  {
    if (!data.active || data.excitation != excitation_idx)
    {
      continue;
    }
    MFEM_VERIFY(std::abs(data.R) > 0.0,
                "Unexpected zero resistance in excited lumped port!");
    for (std::size_t elem_idx = 0; elem_idx < data.elems.size(); elem_idx++)
    {
      const auto &elem = *data.elems[elem_idx];
      const double Rs = data.R * data.GetToSquare(elem);
      const bool terminal_mode = elem_idx < data.terminal_sheet_modes.size() &&
                                 data.terminal_sheet_modes[elem_idx];
      const double Hinc =
          terminal_mode ? std::sqrt(data.R) / Rs
                        : 1.0 / std::sqrt(Rs * elem.GetGeometryWidth() *
                                          elem.GetGeometryLength() * data.elems.size());
      fb.AddCoefficient(data.GetModeCoefficient(elem_idx, 2.0 * Hinc));
    }
  }
}

bool LumpedPortOperator::HasTerminalEdgeExcitation(int excitation_idx) const
{
  (void)excitation_idx;
  return false;
}

void LumpedPortOperator::AddTerminalEdgeExcitationVector(
    int excitation_idx, const mfem::ParFiniteElementSpace &nd_fespace, Vector &rhs) const
{
  Vector lf(nd_fespace.GetVSize());
  lf = 0.0;
  for (const auto &[idx, data] : ports)
  {
    if (!data.active || data.excitation != excitation_idx || !data.HasTerminalEdges())
    {
      continue;
    }
    MFEM_VERIFY(std::abs(data.R) > 0.0,
                "Unexpected zero resistance in excited terminal edge lumped port!");
    data.AddTerminalEdgeExcitationFunctional(nd_fespace, lf, 2.0 / std::sqrt(data.R));
  }
  lf.UseDevice(true);
  nd_fespace.GetProlongationMatrix()->AddMultTranspose(lf, rhs);
}

}  // namespace palace
