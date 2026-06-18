// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_VOLUME_CURRENT_OPERATOR_HPP
#define PALACE_MODELS_VOLUME_CURRENT_OPERATOR_HPP

#include <map>
#include <memory>
#include <vector>
#include <mfem.hpp>

namespace palace
{

class IoData;
class Units;

namespace config
{

struct VolumeCurrentData;

}  // namespace config

//
// Helper class for volume current density sources in a model.
//
class VolumeCurrentData
{
public:
  int excitation;
  mfem::Array<int> attr_list;
  mfem::Vector direction;
  double current_density;
  mfem::Vector source;
  std::unique_ptr<mfem::VectorCoefficient> coef;

public:
  VolumeCurrentData(const config::VolumeCurrentData &data, const mfem::ParMesh &mesh,
                    const Units &units);

  constexpr bool HasExcitation() const { return excitation != 0; }
};

//
// A class handling volume current density sources and their excitation.
//
class VolumeCurrentOperator
{
private:
  // Mapping from source index to data structure containing source current information.
  std::map<int, VolumeCurrentData> sources;

  void SetUpSourceProperties(const std::map<int, config::VolumeCurrentData> &current,
                             const mfem::ParMesh &mesh, const Units &units);
  void PrintSourceInfo(const Units &units);

public:
  VolumeCurrentOperator(const std::map<int, config::VolumeCurrentData> &current,
                        const Units &units, const mfem::ParMesh &mesh);
  VolumeCurrentOperator(const IoData &iodata, const mfem::ParMesh &mesh);

  // Access data structures for the volume current source with the given index.
  const VolumeCurrentData &GetSource(int idx) const;
  auto begin() const { return sources.begin(); }
  auto end() const { return sources.end(); }
  auto rbegin() const { return sources.rbegin(); }
  auto rend() const { return sources.rend(); }
  auto Size() const { return sources.size(); }
  bool Empty() const { return sources.empty(); }

  // Add contributions to the right-hand side source term vector for a volume current
  // excitation, -J_inc for the real version. The frequency-domain iω factor is applied
  // later by SpaceOperator.
  void AddExcitationDomainIntegrator(int idx, mfem::LinearForm &rhs) const;
  void AddExcitationDomainIntegratorsForExcitation(int excitation_idx,
                                                   mfem::LinearForm &rhs) const;
};

}  // namespace palace

#endif  // PALACE_MODELS_VOLUME_CURRENT_OPERATOR_HPP
