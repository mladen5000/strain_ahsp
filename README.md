# Genomic Signature Database

This library provides tools for downloading and analyzing genomic data from NCBI, and for creating a searchable database of genomic signatures.

## Features

- Download genome data from NCBI using their E-utilities API
- Build multi-resolution signatures from genome data
- Store and search signatures in a persistent database
- Support for taxonomic-based searching

## Components

### NCBIDownloader

The `NCBIDownloader` provides functionality for searching and downloading genomes from NCBI:

- `search_genomes`: Search for genomes matching a query
- `download_genome`: Download a genome FASTA file by accession
- `fetch_taxonomy_lineage`: Fetch taxonomic lineage for a taxonomy ID

### SignatureBuilder

The `SignatureBuilder` generates multi-resolution signatures from genome data:

- `build_batch`: Build multiple signatures in batch

### SignatureDatabase

The `SignatureDatabase` provides a persistent database for storing and searching signatures:

- `add_signature`: Add a signature to the database
- `get_signature`: Get a signature by ID
- `search_by_taxonomy`: Search for signatures matching a taxonomy term
- `get_all_signatures`: Get all signatures

### DatabaseManager

The `DatabaseManager` coordinates the downloading, signature generation, and database storage:

- `download_references`: Search for and download reference genomes from NCBI
- `process_references`: Process downloaded genomes into signatures and add to database
- `search_and_add_references`: Search and add reference genomes by query

## Getting Started

### Prerequisites

- Rust 1.56 or later
- Cargo

### Installation

```bash
# Clone the repository
git clone https://github.com/username/genomic-signature-db.git
cd genomic-signature-db

# Build the project
cargo build --release
```

### Usage

```rust
use genomic_signature_db::{DatabaseManager, GenomeMetadata, MultiResolutionSignature};
use std::path::Path;

// Create database manager
let db_path = Path::new("./signature_db");
let cache_dir = Path::new("./genome_cache");
let macro_k = 21;  // k-mer size for macro-resolution
let meso_k = 11;   // k-mer size for meso-resolution
let sketch_size = 1000;  // number of hashes to keep
let threads = 4;  // number of parallel threads
let api_key = None;  // optional NCBI API key

let mut manager = DatabaseManager::new(
    db_path,
    cache_dir,
    macro_k,
    meso_k,
    sketch_size,
    threads,
    api_key,
).unwrap();

// Search and add reference genomes
let added_ids = manager.search_and_add_references("escherichia coli", 5).unwrap();
println!("Added {} reference genomes", added_ids.len());

// Get all signatures
let signatures = manager.database.get_all_signatures().unwrap();
println!("Database contains {} signatures", signatures.len());

// Search by taxonomy
let bacteria_sigs = manager.database.search_by_taxonomy("Bacteria").unwrap();
println!("Found {} bacterial signatures", bacteria_sigs.len());
```

## Testing

The library includes comprehensive unit and integration tests:

```bash
# Run tests
cargo test

# Run tests with output
cargo test -- --nocapture
```

### Unit Tests

Unit tests cover individual components:

- `test_signature_builder_creation`: Test signature builder creation
- `test_signature_builder_batch`: Test batch signature generation
- `test_database_error_conversions`: Test error type conversions
- `test_database_operations`: Test database CRUD operations
- `test_ncbi_downloader`: Test NCBI API interactions
- `test_database_manager`: Test database manager functionality
- `test_genome_metadata_serialization`: Test metadata serialization

### Integration Tests

Integration tests cover multiple components working together:

- `test_download_and_signature_generation`: Test the complete workflow
- `test_signature_database_persistence`: Test database persistence across sessions
- `test_concurrent_database_operations`: Test concurrent database operations
- `test_error_handling`: Test error handling throughout the system

## Dependencies

- `reqwest`: HTTP client for NCBI API access
- `serde`: Serialization/deserialization
- `sled`: Embedded database
- `bincode`: Binary serialization
- `rayon`: Parallel computation
- `log`: Logging
- `thiserror`: Error handling
- `quick-xml`: XML parsing
- `urlencoding`: URL encoding

## Error Handling

The library defines a comprehensive error type `DatabaseError` that wraps errors from various sources:

- IO errors
- HTTP errors
- NCBI API errors
- Database errors
- Serialization errors
- Taxonomy errors
- Signature errors
- Not found errors

## License

This project is licensed under the MIT License - see the LICENSE file for details.
