#pragma once


struct PhysConst
{
    static constexpr double hbarc = 197.3269804;              // MeV·fm
    static constexpr double mN = 939.565420;               // MeV
    static constexpr double two_mc2 = 2.0 * mN;                 // MeV
    static constexpr double hbar2_over_two = (hbarc * hbarc) / two_mc2; // fm²
};