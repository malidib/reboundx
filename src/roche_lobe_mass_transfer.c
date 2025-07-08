/**
 * @file    roche_lobe_mass_transfer.c
 * @brief   Mass transfer between binary stars through Roche-lobe overflow.
 * @author  Generated
 *
 * @section LICENSE
 * Copyright (c) 2024
 *
 * This file is part of reboundx.
 *
 * reboundx is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * reboundx is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 * The section after the dollar signs gets built into the documentation by a script.
 * All lines must start with space * space like below.
 * Tables always must be preceded and followed by a blank line.  See http://docutils.sourceforge.net/docs/user/rst/quickstart.html for a primer on rst.
 * $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 *
 * $Mass Modifications$
 *
 * ======================= ===============================================
 * Author                   Mohamad Ali-Dib
 * Implementation Paper    `Ritter 1988 <https://ui.adsabs.harvard.edu/abs/1988A%26A...202...93R/abstract>`_
 * Based on                `Kolb & Ritter 1990 <https://ui.adsabs.harvard.edu/abs/1990A%26A...236..385K/abstract>`_
 * C Example               :ref:`c_example_roche_lobe_mass_transfer`
 * Python Example          `roche_lobe_mass_transfer.py <https://github.com/dtamayo/reboundx>`_
 * ======================= ===============================================
 *
 * This effect transfers mass from a donor to an accretor star when the donor overfills its Roche lobe.
 * The mass transfer rate follows the prescription of Ritter (1988)
 *
 * .. math::
 *
 *    \dot M = -\dot M_0 \exp\left(\frac{R_d - R_L}{H_p}\right)
 *
 * where :math:`R_d` is the donor radius, :math:`R_L` the Roche-lobe radius (Eggleton 1983) and :math:`H_p`
 * is the photospheric pressure scale height. Mass is removed from the donor and added to the accretor
 * each timestep.
 *
 * **Effect Parameters**
 *
 * ============================ =========== ===============================================
 * Field (C type)               Required    Description
 * ============================ =========== ===============================================
 * rlmt_donor (int)             Yes         Index of donor particle
 * rlmt_accretor (int)          Yes         Index of accretor particle
 * rlmt_loss_fraction (double)  No          Fraction of lost mass that escapes the system
 * ============================ =========== ===============================================
 *
 * **Particle Parameters**
 *
 * ============================ =========== ===============================================
 * Field (C type)               Required    Description
 * ============================ =========== ===============================================
 * rlmt_Hp (double)             Yes         Pressure scale height of the donor star
 * rlmt_mdot0 (double)          Yes         Normalization of the mass transfer rate
 * particles[i].r (double)      Yes         Physical radius of the donor
 * ============================ =========== ===============================================
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"
#include "reboundx.h"

void rebx_roche_lobe_mass_transfer(struct reb_simulation* const sim, struct rebx_operator* const operator, const double dt){
    struct rebx_extras* const rebx = sim->extras;

    int* donor_idx_ptr = rebx_get_param(rebx, operator->ap, "rlmt_donor");
    int* accretor_idx_ptr = rebx_get_param(rebx, operator->ap, "rlmt_accretor");

    if (donor_idx_ptr == NULL || accretor_idx_ptr == NULL){
        rebx_error(rebx, "Need to set rlmt_donor and rlmt_accretor for roche_lobe_mass_transfer.\n");
        return;
    }

    const int donor_idx = *donor_idx_ptr;
    const int accretor_idx = *accretor_idx_ptr;

    if (donor_idx >= sim->N || accretor_idx >= sim->N || donor_idx < 0 || accretor_idx < 0){
        rebx_error(rebx, "rlmt_donor or rlmt_accretor index out of range.\n");
        return;
    }

    struct reb_particle* donor = &sim->particles[donor_idx];
    struct reb_particle* accretor = &sim->particles[accretor_idx];

    const double* Hp_ptr = rebx_get_param(rebx, donor->ap, "rlmt_Hp");
    const double* mdot0_ptr = rebx_get_param(rebx, donor->ap, "rlmt_mdot0");
    const double* loss_frac_ptr = rebx_get_param(rebx, operator->ap, "rlmt_loss_fraction");

    if (Hp_ptr == NULL || mdot0_ptr == NULL){
        rebx_error(rebx, "Need to set rlmt_Hp and rlmt_mdot0 on donor particle.\n");
        return;
    }

    const double Hp = *Hp_ptr;
    const double mdot0 = *mdot0_ptr;
    double loss_frac = 0.;
    if (loss_frac_ptr != NULL){
        loss_frac = *loss_frac_ptr;
    }
    if (loss_frac < 0.){
        loss_frac = 0.;
    } else if (loss_frac > 1.){
        loss_frac = 1.;
    }

    const double dx = donor->x - accretor->x;
    const double dy = donor->y - accretor->y;
    const double dz = donor->z - accretor->z;
    const double r = sqrt(dx*dx + dy*dy + dz*dz);

    const double q = donor->m / accretor->m;
    const double q13 = cbrt(q);
    const double rl = r * (0.49*q13*q13) / (0.6*q13*q13 + log(1.+q13));
    const double Rd = donor->r;


    const double mdot = -mdot0 * exp((Rd - rl)/Hp); // negative for mass loss
    double dM = mdot * dt;
    if (donor->m + dM <= 0.){
        dM = -donor->m; // remove remaining mass only once
    }

    const double mass_loss = -dM; // positive amount of mass leaving donor
    const double mass_accreted = mass_loss * (1. - loss_frac);

    donor->m -= mass_loss;
    accretor->m += mass_accreted;

    // Stop accretion once the donor has lost all of its mass. The donor
    // particle is removed from the simulation and the operator itself is
    // detached so it no longer runs.

    if (donor->m <= 0.){
        reb_simulation_remove_particle(sim, donor_idx, 1);
        rebx_remove_operator(rebx, operator);
        reb_simulation_move_to_com(sim);
        return;
    }


    reb_simulation_move_to_com(sim);
}

