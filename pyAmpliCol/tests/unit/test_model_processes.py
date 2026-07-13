from __future__ import annotations

import pytest

from pyamplicol.model_assets import bundled_model_path
from pyamplicol.model_processes import (
    ModelParticleCatalog,
    build_model_process_ir,
    expand_model_processes,
    parse_multiparticle_definitions,
)
from pyamplicol.model_source import compile_model_source


def test_builtin_catalog_contains_particle_and_antiparticle_records() -> None:
    compiled = compile_model_source("BUILTIN_SM", use_cache=False)
    catalog = ModelParticleCatalog(compiled.name, compiled.ir.particles)

    process = build_model_process_ir("d d~ > Z g", compiled.ir)

    assert catalog.resolve("d").antiname == "d~"
    assert catalog.resolve("d~").antiname == "d"
    assert process.process == "d d~ > z g"
    assert process.outgoing_pdgs == (-1, 1, 23, 21)
    assert process.quark_lines.quark_pair_count == 1


def test_ufo_catalog_resolves_case_only_when_unambiguous() -> None:
    compiled = compile_model_source(
        bundled_model_path("sm", "json"),
        use_cache=False,
    )

    process = build_model_process_ir("d d~ > z h", compiled.ir)

    assert process.process == "d d~ > Z H"
    assert process.final_pdgs == (23, 25)


def test_only_sm_catalogs_ship_p_and_j_aliases() -> None:
    sm = compile_model_source(
        bundled_model_path("sm", "json"),
        use_cache=False,
    )
    scalars = compile_model_source(
        bundled_model_path("scalars", "json"),
        use_cache=False,
    )

    sm_aliases = ModelParticleCatalog(sm.name, sm.ir.particles).default_multiparticles()
    scalar_aliases = ModelParticleCatalog(
        scalars.name,
        scalars.ir.particles,
    ).default_multiparticles()

    assert set(sm_aliases) == {"p", "j"}
    assert "g" in sm_aliases["p"]
    assert scalar_aliases == {}


def test_user_multiparticle_definitions_expand_repeated_slots() -> None:
    compiled = compile_model_source(
        bundled_model_path("scalars", "json"),
        use_cache=False,
    )
    catalog = ModelParticleCatalog(compiled.name, compiled.ir.particles)
    aliases = parse_multiparticle_definitions(
        ["phi=scalar_0,scalar_1"],
        catalog,
    )

    expanded = expand_model_processes(
        "scalar_0 scalar_0 > 2*phi",
        catalog,
        multiparticles=aliases,
    )

    assert expanded == (
        "scalar_0 scalar_0 > scalar_0 scalar_0",
        "scalar_0 scalar_0 > scalar_0 scalar_1",
        "scalar_0 scalar_0 > scalar_1 scalar_0",
        "scalar_0 scalar_0 > scalar_1 scalar_1",
    )


def test_catalog_rejects_ghosts_and_unknown_particles_as_external_states() -> None:
    compiled = compile_model_source(
        bundled_model_path("sm", "json"),
        use_cache=False,
    )
    catalog = ModelParticleCatalog(compiled.name, compiled.ir.particles)

    with pytest.raises(ValueError, match="not an external state"):
        catalog.resolve("ghG")
    with pytest.raises(ValueError, match="not an external state"):
        catalog.resolve("not-a-particle")
