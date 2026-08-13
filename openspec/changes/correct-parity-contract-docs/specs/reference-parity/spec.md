# Reference Parity Change Specification

## ADDED Requirements

### Requirement: Parity claims identify reference prerequisites
Documentation SHALL distinguish Rust-only tests from tests that require a C minibwa binary or external reference data.

#### Scenario: User reads validation instructions
- **WHEN** the user reads the test and compatibility sections
- **THEN** each live parity command SHALL list its required binary and data
- **AND** SHALL state that missing prerequisites mean parity was not verified

### Requirement: Current CI scope is truthful
Documentation SHALL NOT claim that current CI proves C parity while its parity tests can skip for missing prerequisites.

#### Scenario: CI runs without C reference
- **GIVEN** the workflow does not build/provide C minibwa and required data
- **WHEN** Cargo tests complete
- **THEN** project documentation SHALL classify the run as Rust test success, not C parity success

### Requirement: Skip behavior is observable
When a live parity test is skipped because a reference dependency is absent, the reason SHALL be visible in test output or documented verification evidence.

#### Scenario: Reference binary missing
- **WHEN** a developer runs the parity suite without the configured C binary
- **THEN** the run SHALL expose that live parity was not executed
- **AND** SHALL NOT describe the skipped case as a passed comparison

### Requirement: Comparison helper matches supported interfaces
Any maintained comparison script SHALL invoke commands and consume files that exist in the audited Rust and C interfaces.

#### Scenario: Helper cannot be validated
- **WHEN** the C source/version needed to validate the helper is unknown
- **THEN** the helper change SHALL remain blocked or the script SHALL be explicitly documented as unsupported
- **AND** no substitute implementation SHALL be presented as the reference
