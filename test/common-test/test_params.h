#pragma once

constexpr const float float_test_epsilon = 1e-5;
// reprojection error can be more than epsilon
const constexpr float reprojection_tol = 5.f;
const constexpr float svd_tol = 0.01f;
const constexpr float localization_tol_rel = 0.05f;
const constexpr float localization_tol_abs = 0.1f;
