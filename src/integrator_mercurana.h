/**
 * @file 	integrator_mercurana.h
 * @brief 	Interface for numerical particle integrator
 * @author 	Hanno Rein 
 * 
 * @section 	LICENSE
 * Copyright (c) 2017 Hanno Rein
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
#ifndef _INTEGRATOR_MERCURANA_H
#define _INTEGRATOR_MERCURANA_H
void reb_integrator_mercurana_part1(struct reb_simulation* r);          ///< Internal function used to call a specific integrator
void reb_integrator_mercurana_part2(struct reb_simulation* r);          ///< Internal function used to call a specific integrator
void reb_integrator_mercurana_synchronize(struct reb_simulation* r);    ///< Internal function used to call a specific integrator
void reb_integrator_mercurana_reset(struct reb_simulation* r);          ///< Internal function used to call a specific integrator
void reb_integrator_mercurana_inertial_to_dh(struct reb_simulation* r); ///< Internal in-place coordinate transformation
#endif
