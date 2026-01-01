# Leinster Ring Counter

A Python implementation for counting and analyzing Leinster rings, based on the mathematical theory from "The Arithmetic of Ideal Lattices".

## Overview

Leinster rings are a special class of finite rings that satisfy certain algebraic conditions related to their ideal structure. This project provides tools to enumerate, classify, and analyze these rings up to given orders.

## Features

- **Ring Generation**: Generate all finite rings of given orders
- **Leinster Detection**: Identify which rings satisfy the Leinster condition
- **Classification**: Categorize rings into Type I, II, and III
- **Parallel Processing**: Multi-core optimization for large-scale counting
- **Pattern Recognition**: Analyze structural patterns in Leinster rings
- **Advanced Pruning**: Efficient elimination of impossible orders
- **Reporting**: Generate text and CSV reports of results

## Installation

1. Clone the repository:
```bash
git clone https://github.com/Sulaiman-Sabo/Linster-Ring-Project.git
cd Linster-Ring-Project
```

2. Create a virtual environment:
```bash
python -m venv .venv
.venv\Scripts\activate  # On Windows
```

3. Install dependencies:
```bash
pip install numpy sympy
```

## Usage

### Basic Counting

```python
from LinsterRingCounter import main

# Count Leinster rings up to order 100
result, analysis = main()
```

### Advanced Features

```python
from LinsterRingCounter import ParallelCounter, AdvancedPruning, PatternRecognizer

# Parallel counting
counter = ParallelCounter(max_order=1000, num_processes=4)
result = counter.count_parallel()

# Advanced pruning
pruner = AdvancedPruning()
possible_orders = pruner.generate_possible_orders(1000)

# Pattern analysis
analyzer = PatternRecognizer()
patterns = analyzer.analyze_patterns(result["rings"])
```

## Project Structure

```
Linster-Ring-Project/
├── LinsterRingCounter.py    # Main implementation
├── backup.py                # Backup of previous version
├── leinster_rings_10.csv    # Sample output data
├── docs/                    # Documentation
│   └── description.md       # Detailed project description
├── README.md                # This file
└── .venv/                   # Virtual environment (not in repo)
```

## Mathematical Background

Leinster rings are defined by the condition that the sum of the reciprocals of the sizes of their ideals equals 2. This project implements algorithms to:

1. Generate finite rings using group theory and ring constructions
2. Compute the σ-function (sum of ideal sizes)
3. Check the Leinster condition: σ(R) = 2|R|
4. Classify rings based on their construction

## Contributing

Contributions are welcome! Please feel free to submit issues, feature requests, or pull requests.

## License

This project is open source. Please check the license file for details.

## References

- Leinster, T. "The Arithmetic of Ideal Lattices" (research paper)
- Ring theory and algebraic structures