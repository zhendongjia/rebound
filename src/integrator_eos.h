/**
 * @file 	integrator_eos.h
 * @brief 	Interface for numerical particle integrator
 * @author 	Hanno Rein 
 * 
 * @section 	LICENSE
 * Copyright (c) 2019 Hanno Rein
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#ifndef _INTEGRATOR_EOS_H
#define _INTEGRATOR_EOS_H
void reb_integrator_eos_part1(struct reb_simulation* r);          ///< Internal function used to call a specific integrator
void reb_integrator_eos_part2(struct reb_simulation* r);          ///< Internal function used to call a specific integrator
void reb_integrator_eos_synchronize(struct reb_simulation* r);    ///< Internal function used to call a specific integrator
void reb_integrator_eos_reset(struct reb_simulation* r);          ///< Internal function used to call a specific integrator

void reb_integrator_eos_step(struct reb_simulation* const r, double dt, double dtfacfirst, double dtfaclast, int shell, enum REB_EOS_TYPE type, void (*drift_step)(struct reb_simulation* const r, double a, unsigned int shell), void (*interaction_step)(struct reb_simulation* const r, double y, double v, unsigned int shell));
void reb_integrator_eos_preprocessor(struct reb_simulation* const r, double dt, unsigned int shell, enum REB_EOS_TYPE type, void (*drift_step)(struct reb_simulation* const r, double a, unsigned int shell), void (*interaction_step)(struct reb_simulation* const r, double y, double v, unsigned int shell));
void reb_integrator_eos_postprocessor(struct reb_simulation* const r, double dt, unsigned int shell, enum REB_EOS_TYPE type, void (*drift_step)(struct reb_simulation* const r, double a, unsigned int shell), void (*interaction_step)(struct reb_simulation* const r, double y, double v, unsigned int shell));

extern const double reb_eos_lf4_a;
extern const double reb_eos_lf6_a[5];
extern const double reb_eos_lf8_a[9];
extern const double reb_eos_lf4_2_a;
extern const double reb_eos_lf8_6_4_a[4];
extern const double reb_eos_lf8_6_4_b[4];
extern const double reb_eos_pmlf6_a[2];
extern const double reb_eos_pmlf6_b[2];
extern const double reb_eos_pmlf6_c[2];
extern const double reb_eos_pmlf6_z[6];
extern const double reb_eos_pmlf6_y[6];
extern const double reb_eos_pmlf6_v[6];
extern const double reb_eos_pmlf4_y[3];
extern const double reb_eos_pmlf4_z[3];
extern const double reb_eos_plf7_6_4_a[2];
extern const double reb_eos_plf7_6_4_b[2];
extern const double reb_eos_plf7_6_4_z[6];
extern const double reb_eos_plf7_6_4_y[6];

#endif
