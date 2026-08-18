# Reference Parity Change Specification

## ADDED Requirements

### Requirement: Parity claims identify reference prerequisites
Documentation SHALL distinguish Rust-only tests from tests that require a C minibwa binary or external reference data.

#### Scenario: User reads validation instructions
- **WHEN** the user reads the test and compatibility sections
- **THEN** each live parity command SHALL list its required binary and data
- **AND** SHALL state that missing prerequisites mean parity was not verified

### Requirement: CI builds the pinned C reference
CI SHALL clone `https://github.com/lh3/minibwa` at the commit in `C_MINIBWA_COMMIT`, build `minibwa`, place it at `../minibwa/minibwa` relative to the crate root, and set `MINIBWA_REQUIRE_PARITY=1`.

#### Scenario: CI runs with C reference
- **GIVEN** the workflow builds the pinned C minibwa and test data from that clone
- **WHEN** Cargo tests complete
- **THEN** live parity tests SHALL have executed
- **AND** a missing binary or chrM fixture SHALL fail the job

### Requirement: Local skip remains observable
When a developer runs tests without the C binary and without `MINIBWA_REQUIRE_PARITY`, skip SHALL print `[PARITY SKIP]` and SHALL NOT be described as a passed comparison.

#### Scenario: Reference binary missing locally
- **WHEN** a developer runs the parity suite without `../minibwa/minibwa` and without `MINIBWA_REQUIRE_PARITY`
- **THEN** the run SHALL expose that live parity was not executed
- **AND** SHALL NOT describe the skipped case as a passed comparison

### Requirement: Comparison helper matches supported interfaces
`scripts/compare.sh` SHALL invoke C `minibwa map` and Rust `map` against `.l2b`/`.mbw` indexes produced by each binary's `index` command.

#### Scenario: Helper is run with current CLI
- **WHEN** the user runs `./scripts/compare.sh ref.fa reads.fq`
- **THEN** the script SHALL fail if the C binary is missing
- **AND** SHALL NOT invoke `bwa` or a non-existent Rust `align` subcommand
