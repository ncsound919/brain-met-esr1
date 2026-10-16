"""Two-state stochastic dormancy model: dormant (D) <-> active (A).

Bet-hedging prediction: pause-targeting the active state increases
diversification risk. GillesPy2 API: rates are Parameter objects.
"""

from __future__ import annotations


def build_model(
    switching_rate: float = 0.01,
    reverse_rate: float | None = None,
    proliferation_rate: float = 0.1,
    death_rate: float = 0.05,
    pause_target_rate: float = 0.0,
    d0: int = 100,
    a0: int = 0,
):
    import gillespy2

    reverse_rate = switching_rate if reverse_rate is None else reverse_rate
    model = gillespy2.Model(name="dormancy_bet_hedging")

    k_switch = gillespy2.Parameter(name="k_switch", expression=switching_rate)
    k_rev = gillespy2.Parameter(name="k_rev", expression=reverse_rate)
    k_prolif = gillespy2.Parameter(name="k_prolif", expression=proliferation_rate)
    k_death = gillespy2.Parameter(name="k_death", expression=death_rate)
    model.add_parameter([k_switch, k_rev, k_prolif, k_death])

    D = gillespy2.Species(name="D", initial_value=d0)
    A = gillespy2.Species(name="A", initial_value=a0)
    model.add_species([D, A])

    model.add_reaction(
        gillespy2.Reaction(name="switch", reactants={D: 1}, products={A: 1}, rate=k_switch)
    )
    model.add_reaction(
        gillespy2.Reaction(name="reverse", reactants={A: 1}, products={D: 1}, rate=k_rev)
    )
    model.add_reaction(
        gillespy2.Reaction(name="proliferate", reactants={A: 1}, products={A: 2}, rate=k_prolif)
    )
    model.add_reaction(
        gillespy2.Reaction(name="death", reactants={A: 1}, products={}, rate=k_death)
    )

    if pause_target_rate > 0:
        k_pause = gillespy2.Parameter(name="k_pause", expression=pause_target_rate)
        model.add_parameter([k_pause])
        model.add_reaction(
            gillespy2.Reaction(name="pause_target", reactants={A: 1}, products={}, rate=k_pause)
        )
    model.timespan(range(501))
    return model
