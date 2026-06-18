// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "volumecurrentoperator.hpp"

#include <algorithm>
#include <fmt/ranges.h>
#include "fem/coefficient.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

namespace palace
{

namespace
{

double CurrentDensityScale(const Units &units)
{
  const double current_scale = units.GetScaleFactor<Units::ValueType::CURRENT>();
  const double length_scale = units.GetScaleFactor<Units::ValueType::LENGTH>();
  return current_scale / (length_scale * length_scale);
}

}  // namespace

VolumeCurrentData::VolumeCurrentData(const config::VolumeCurrentData &data,
                                     const mfem::ParMesh &mesh, const Units &units)
  : excitation(data.excitation), direction(mesh.SpaceDimension()),
    current_density(data.current_density / CurrentDensityScale(units)),
    source(mesh.SpaceDimension())
{
  attr_list.Append(data.attributes.data(), data.attributes.size());

  MFEM_VERIFY(data.direction.size() == static_cast<std::size_t>(mesh.SpaceDimension()),
              "Volume current direction dimension must match mesh space dimension!");
  for (int d = 0; d < mesh.SpaceDimension(); d++)
  {
    direction[d] = data.direction[d];
  }

  source = direction;
  source *= -current_density;
  coef = std::make_unique<RestrictedVectorCoefficient<mfem::VectorConstantCoefficient>>(
      attr_list, source);
}

VolumeCurrentOperator::VolumeCurrentOperator(const IoData &iodata,
                                             const mfem::ParMesh &mesh)
{
  SetUpSourceProperties(iodata.domains.volume_current, mesh, iodata.units);
  PrintSourceInfo(iodata.units);
}

VolumeCurrentOperator::VolumeCurrentOperator(
    const std::map<int, config::VolumeCurrentData> &current, const Units &units,
    const mfem::ParMesh &mesh)
{
  SetUpSourceProperties(current, mesh, units);
  PrintSourceInfo(units);
}

void VolumeCurrentOperator::SetUpSourceProperties(
    const std::map<int, config::VolumeCurrentData> &current, const mfem::ParMesh &mesh,
    const Units &units)
{
  if (!current.empty())
  {
    int attr_max = mesh.attributes.Size() ? mesh.attributes.Max() : 0;
    mfem::Array<int> attr_marker(attr_max), source_marker(attr_max);
    attr_marker = 0;
    source_marker = 0;
    for (auto attr : mesh.attributes)
    {
      attr_marker[attr - 1] = 1;
    }
    for (const auto &[idx, data] : current)
    {
      for (auto attr : data.attributes)
      {
        MFEM_VERIFY(attr > 0 && attr <= attr_max,
                    "Volume current source attribute tags must be non-negative and "
                    "correspond to domain attributes in the mesh!");
        MFEM_VERIFY(attr_marker[attr - 1],
                    "Unknown volume current source attribute " << attr << "!");
        MFEM_VERIFY(
            !source_marker[attr - 1],
            "Domain attribute is assigned to more than one volume current source!");
        source_marker[attr - 1] = 1;
      }
    }
  }

  for (const auto &[idx, data] : current)
  {
    sources.try_emplace(idx, data, mesh, units);
  }
}

void VolumeCurrentOperator::PrintSourceInfo(const Units &units)
{
  if (sources.empty())
  {
    return;
  }

  Mpi::Print("\nConfiguring volume current excitation source terms:\n");
  for (const auto &[idx, data] : sources)
  {
    const double physical_current_density =
        data.current_density * CurrentDensityScale(units);
    Mpi::Print(" Source {:d}: Excitation = {:d}, Attributes = [{}]\n"
               " \tCurrentDensity = {:.3e} A/m^2\n"
               " \tDirection = ({:.3e})\n",
               idx, data.excitation, fmt::join(data.attr_list, ", "),
               physical_current_density, fmt::join(data.direction, ", "));
  }
}

const VolumeCurrentData &VolumeCurrentOperator::GetSource(int idx) const
{
  auto it = sources.find(idx);
  MFEM_VERIFY(it != sources.end(), "Unknown volume current source index requested!");
  return it->second;
}

void VolumeCurrentOperator::AddExcitationDomainIntegrator(int idx,
                                                          mfem::LinearForm &rhs) const
{
  const auto &source = GetSource(idx);
  rhs.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(*source.coef));
}

void VolumeCurrentOperator::AddExcitationDomainIntegratorsForExcitation(
    int excitation_idx, mfem::LinearForm &rhs) const
{
  for (const auto &[idx, source] : sources)
  {
    if (source.excitation == excitation_idx)
    {
      rhs.AddDomainIntegrator(new mfem::VectorFEDomainLFIntegrator(*source.coef));
    }
  }
}

}  // namespace palace
