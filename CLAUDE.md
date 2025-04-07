# Metagenomics DESeq2 Development Guide

## Build & Test Commands
```bash
cargo build                       # Build in debug mode
cargo build --release             # Build in release mode
cargo build --example <example>   # Build specific example
cargo test                        # Run all tests
cargo test <test_name>            # Run specific test
cargo clippy -- -D warnings       # Run linter
cargo fmt                         # Format code
cargo fmt -- --check              # Check format without changes
```

## Code Style Guidelines
- **Naming**: `snake_case` for variables/functions/modules, `CamelCase` for types/traits
- **Imports**: Group by 1) std library 2) external crates 3) local imports
- **Error Handling**: Use `anyhow::Result` for functions with multiple error types and `thiserror::Error` for custom errors
- **Types**: Prefer strong typing with custom structs over primitives; use `Option<T>` and `Result<T, E>` appropriately
- **Documentation**: Document all public items with `///` comments; include description, arguments, returns
- **Tests**: Place in `#[cfg(test)]` module; use descriptive names with `test_` prefix
- **Organization**: Group by functionality; follow single responsibility principle

This project uses Rust 2021 edition with standard formatting and modules organized by domain functionality.