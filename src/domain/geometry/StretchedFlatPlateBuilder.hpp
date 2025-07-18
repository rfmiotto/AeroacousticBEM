#pragma once

#include "FlatPlate.hpp"
#include "Stretching.hpp"

namespace bem::domain::geometry {

FlatPlate2D makeStretchedPlate(std::size_t nAirfoilElements,
                               std::size_t nVerticalElements,
                               Real thickness,
                               Real domainLength,
                               const StretchingParams &params);

} // namespace bem::domain::geometry
