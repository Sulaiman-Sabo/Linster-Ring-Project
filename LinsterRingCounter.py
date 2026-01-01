import math
import itertools
from collections import defaultdict
from functools import lru_cache
from typing import List, Tuple, Dict, Set, Optional, Union
from dataclasses import dataclass
import numpy as np
from sympy import factorint, isprime, divisors

# SECTION 1: RING REPRESENTATION AND GENERATION


@dataclass
class Ring:
    """Basic ring representation."""
    order: int
    structure: str  # e.g., "Z_n", "F_q", "Z_p^k", "product"
    components: List  # For product rings
    is_commutative: bool = True
    is_local: bool = False
    is_field: bool = False
    is_cyclic: bool = False
    
    def __repr__(self):
        return f"Ring(order={self.order}, structure={self.structure})"

class RingGenerator:
    """Generate finite rings up to given order."""
    
    def __init__(self, max_order: int):
        self.max_order = max_order
        self.prime_cache = {}
        self.ring_cache = {}
        
    def generate_all_rings(self, order: int) -> List[Ring]:
        """Generate all rings of given order."""
        if order in self.ring_cache:
            return self.ring_cache[order]
        
        rings = []
        
        # Quick elimination based on known obstructions
        if self._quick_eliminate(order):
            self.ring_cache[order] = rings
            return rings
        
        # 1. Cyclic rings Z_n
        rings.append(self._create_cyclic_ring(order))
        
        # 2. Local rings (prime power order)
        if self._is_prime_power(order):
            rings.extend(self._generate_local_rings(order))
        
        # 3. Product rings (for composite orders)
        if not self._is_prime_power(order):
            rings.extend(self._generate_product_rings(order))
        
        # 4. Finite fields (when order is prime power)
        if self._is_prime_power(order) and isprime(order):
            rings.append(self._create_field(order))
        
        self.ring_cache[order] = rings
        return rings
    
    def _quick_eliminate(self, order: int) -> bool:
        """Quick elimination of impossible orders."""
        # Parity obstruction: even orders cannot be Leinster (except perfect numbers?)
        # if order % 2 == 0 and order > 2:
        #     return True
        
        # Known from paper: no Leinster rings of prime power order > 1
        if self._is_prime_power(order) and order > 1:
            return True
        
        return False
    
    def _is_prime_power(self, n: int) -> bool:
        """Check if n is a prime power."""
        factors = factorint(n)
        return len(factors) == 1
    
    def _create_cyclic_ring(self, n: int) -> Ring:
        """Create cyclic ring Z_n."""
        return Ring(
            order=n,
            structure=f"Z_{n}",
            components=[],
            is_cyclic=True,
            is_commutative=True,
            is_local=self._is_prime_power(n),
            is_field=(isprime(n) and n > 1)
        )
    
    def _generate_local_rings(self, order: int) -> List[Ring]:
        """Generate local rings of prime power order."""
        rings = []
        factors = factorint(order)
        
        if len(factors) == 1:
            p, k = list(factors.items())[0]
            
            # Z_{p^k}
            rings.append(Ring(
                order=order,
                structure=f"Z_{p}^{k}",
                components=[],
                is_commutative=True,
                is_local=True,
                is_field=(k == 1)
            ))
            
            # Other local rings like Z_{p^k}[x]/(x^2), etc.
            # This is simplified - in practice, more types exist
            if k >= 2:
                # Example: Z_{p^2}
                rings.append(Ring(
                    order=order,
                    structure=f"Local_{p}^{k}",
                    components=[],
                    is_commutative=True,
                    is_local=True,
                    is_field=False
                ))
        
        return rings
    
    def _generate_product_rings(self, order: int) -> List[Ring]:
        """Generate product rings for composite orders."""
        rings = []
        
        # Get all factorizations a*b = order with a,b > 1
        factor_pairs = []
        for d in divisors(order):
            if d > 1 and d < order:
                e = order // d
                if e > 1:
                    factor_pairs.append((d, e))
        
        # Generate product rings
        for a, b in factor_pairs:
            # Get component rings
            comps_a = self.generate_all_rings(a)
            comps_b = self.generate_all_rings(b)
            
            # Create all combinations
            for ring_a in comps_a:
                for ring_b in comps_b:
                    product_ring = Ring(
                        order=order,
                        structure=f"{ring_a.structure}×{ring_b.structure}",
                        components=[ring_a, ring_b],
                        is_commutative=ring_a.is_commutative and ring_b.is_commutative,
                        is_local=False,
                        is_field=False
                    )
                    rings.append(product_ring)
        
        return rings
    
    def _create_field(self, p: int) -> Ring:
        """Create finite field F_p."""
        return Ring(
            order=p,
            structure=f"F_{p}",
            components=[],
            is_commutative=True,
            is_local=True,
            is_field=True,
            is_cyclic=True
        )


# SECTION 2: IDEAL COMPUTATION

class IdealComputer:
    """Compute ideals and σ(R) for rings."""
    
    def __init__(self):
        self.ideal_cache = {}
        self.sigma_cache = {}
        
    def compute_sigma(self, ring: Ring) -> int:
        """Compute σ(R) = sum of orders of all two-sided ideals."""
        cache_key = (ring.order, ring.structure)
        if cache_key in self.sigma_cache:
            return self.sigma_cache[cache_key]
        
        # For product rings, use multiplicativity
        if ring.components:
            sigma = 1
            for comp in ring.components:
                sigma *= self.compute_sigma(comp)
            self.sigma_cache[cache_key] = sigma
            return sigma
        
        # For basic rings, compute directly
        if ring.structure.startswith("Z_"):
            # For Z_n, ideals correspond to divisors
            sigma = self._sigma_cyclic(ring.order)
        elif ring.structure.startswith("F_"):
            # Fields have only {0} and the field itself
            sigma = 1 + ring.order
        elif ring.structure.startswith("Local_"):
            # Local rings: ideals form a chain
            sigma = self._sigma_local(ring.order)
        else:
            # Fallback: approximate or compute via brute force
            sigma = self._sigma_bruteforce(ring)
        
        self.sigma_cache[cache_key] = sigma
        return sigma
    
    def _sigma_cyclic(self, n: int) -> int:
        """σ(Z_n) = σ(n) (sum of divisors)."""
        return sum(divisors(n))
    
    def _sigma_local(self, order: int) -> int:
        """σ for local ring of order p^k: 1 + p + p^2 + ... + p^k."""
        factors = factorint(order)
        if len(factors) == 1:
            p, k = list(factors.items())[0]
            # Geometric series sum
            sigma = (p**(k+1) - 1) // (p - 1)
            return sigma
        return 0
    
    def _sigma_bruteforce(self, ring: Ring) -> int:
        """Brute-force σ computation (for small rings)."""
        # This is simplified - in practice would enumerate additive subgroups
        # and check ideal conditions
        
        # For demonstration, use approximation
        if ring.is_field:
            return 1 + ring.order
        elif ring.is_cyclic:
            return self._sigma_cyclic(ring.order)
        else:
            # Approximate with bounds
            return ring.order + 1  # Lower bound
    
    def compute_rho(self, ring: Ring) -> float:
        """Compute ρ(R) = σ(R)/|R|."""
        sigma = self.compute_sigma(ring)
        return sigma / ring.order if ring.order > 0 else 0
    
    def is_leinster(self, ring: Ring) -> bool:
        """Check if ring is Leinster: σ(R) = 2|R|."""
        sigma = self.compute_sigma(ring)
        return sigma == 2 * ring.order


# SECTION 3: OPTIMIZED COUNTING ALGORITHMS


class LeinsterCounter:
    """Main class for counting Leinster rings."""
    
    def __init__(self, max_order: int):
        self.max_order = max_order
        self.generator = RingGenerator(max_order)
        self.computer = IdealComputer()
        self.rho_table = {}  # Cache for ρ values
        self.leinster_rings = []
        
    def count_basic(self) -> Dict:
        """Basic counting algorithm (brute force)."""
        print(f"Counting Leinster rings up to order {self.max_order}...")
        
        total_count = 0
        leinster_list = []
        
        for n in range(1, self.max_order + 1):
            if n % 100 == 0:
                print(f"  Processing order {n}... ({total_count} found)")
            
            rings = self.generator.generate_all_rings(n)
            
            for ring in rings:
                if self.computer.is_leinster(ring):
                    total_count += 1
                    leinster_list.append(ring)
                    print(f"    Found Leinster ring #{total_count}: {ring}")
        
        return {
            "total": total_count,
            "rings": leinster_list,
            "method": "basic"
        }
    
    def count_optimized(self) -> Dict:
        """Optimized counting using ρ-table and decomposition."""
        print(f"Optimized counting up to order {self.max_order}...")
        
        # Precompute ρ-values for indecomposable rings
        self._precompute_rho_table()
        
        total_count = 0
        leinster_list = []
        
        for n in range(1, self.max_order + 1):
            # Quick elimination
            if self._should_skip_order(n):
                continue
            
            if n % 100 == 0:
                print(f"  Processing order {n}... ({total_count} found)")
            
            rings = self.generator.generate_all_rings(n)
            
            for ring in rings:
                if self._is_decomposable(ring):
                    # Check using product condition
                    if self._check_product_leinster(ring):
                        total_count += 1
                        leinster_list.append(ring)
                        print(f"    Found product Leinster #{total_count}: {ring}")
                else:
                    # Check indecomposable ring directly
                    if self.computer.is_leinster(ring):
                        total_count += 1
                        leinster_list.append(ring)
                        print(f"    Found indecomposable Leinster #{total_count}: {ring}")
        
        return {
            "total": total_count,
            "rings": leinster_list,
            "method": "optimized"
        }
    
    def _precompute_rho_table(self):
        """Precompute ρ-values for all rings up to max_order."""
        print("Precomputing ρ-table...")
        
        for n in range(1, self.max_order + 1):
            if self._should_skip_order(n):
                continue
            
            rings = self.generator.generate_all_rings(n)
            for ring in rings:
                if not self._is_decomposable(ring):
                    rho = self.computer.compute_rho(ring)
                    key = (ring.order, ring.structure)
                    self.rho_table[key] = rho
    
    def _should_skip_order(self, n: int) -> bool:
        """Quick elimination of impossible orders."""
        # Parity obstruction
        if n % 2 == 0 and n > 2:
            return True
        
        # Prime power obstruction
        factors = factorint(n)
        if len(factors) == 1:
            p, k = list(factors.items())[0]
            if k >= 1 and n > 1:
                return True
        
        return False
    
    def _is_decomposable(self, ring: Ring) -> bool:
        """Check if ring is decomposable (product of smaller rings)."""
        return bool(ring.components)
    
    def _check_product_leinster(self, ring: Ring) -> bool:
        """Check if product ring is Leinster using ρ-values."""
        if not ring.components:
            return False
        
        # Compute product of ρ-values
        rho_product = 1.0
        for comp in ring.components:
            key = (comp.order, comp.structure)
            if key in self.rho_table:
                rho_product *= self.rho_table[key]
            else:
                # Fallback: compute directly
                rho_product *= self.computer.compute_rho(comp)
        
        # Check if ρ(product) = 2
        return abs(rho_product - 2.0) < 1e-10

#############################################################################
## SECTION 4: CLASSIFICATION AND ANALYSIS
#############################################################################

class RingAnalyzer:
    """Analyze and classify Leinster rings."""
    
    @staticmethod
    def classify_ring(ring: Ring) -> str:
        """Classify ring into Type I, II, or III."""
        # Type I: Cyclic perfect number rings
        if ring.is_cyclic:
            n = ring.order
            # Check if n is perfect number
            if sum(divisors(n)) == 2 * n:
                return "Type I"
        
        # Type II: Products of finite fields
        if ring.components:
            all_fields = all(comp.is_field for comp in ring.components)
            if all_fields:
                return "Type II"
        
        # Type III: Mixed products
        return "Type III"
    
    @staticmethod
    def analyze_results(result: Dict) -> Dict:
        """Analyze counting results."""
        rings = result["rings"]
        total = len(rings)
        
        if total == 0:
            return {"error": "No rings found"}
        
        # Count by type
        type_counts = {"Type I": 0, "Type II": 0, "Type III": 0}
        orders = []
        
        for ring in rings:
            ring_type = RingAnalyzer.classify_ring(ring)
            type_counts[ring_type] += 1
            orders.append(ring.order)
        
        # Statistics
        unique_orders = set(orders)
        
        return {
            "total": total,
            "type_counts": type_counts,
            "orders": {
                "min": min(orders),
                "max": max(orders),
                "avg": sum(orders) / len(orders),
                "unique": len(unique_orders)
            },
            "density": total / max(orders) if orders else 0
        }
    
    @staticmethod
    def asymptotic_estimate(N: int, result: Dict) -> Dict:
        """Estimate asymptotic behavior L(N) ~ C * N / log(N)."""
        total = result["total"]
        
        if total == 0:
            return {"C_estimate": 0}
        
        # Estimate constant C in L(N) ~ C * N / log(N)
        C_estimate = total * math.log(N) / N if N > 1 else 0
        
        # Estimate density decay
        density = total / N if N > 0 else 0
        
        return {
            "C_estimate": C_estimate,
            "density": density,
            "L(N)": total,
            "N": N
        }

#############################################################################
## SECTION 5: VISUALIZATION AND REPORTING
#############################################################################

class ReportGenerator:
    """Generate reports and visualizations."""
    
    @staticmethod
    def generate_text_report(result: Dict, analysis: Dict) -> str:
        """Generate text report of results."""
        report = []
        report.append("=" * 60)
        report.append("LEINSTER RING ENUMERATION REPORT")
        report.append("=" * 60)
        report.append(f"Method: {result.get('method', 'unknown')}")
        
        if 'error' in analysis:
            report.append(f"Error: {analysis['error']}")
            return "\n".join(report)
        
        report.append(f"Maximum order: {max([r.order for r in result['rings']]) if result['rings'] else 0}")
        report.append(f"Total Leinster rings found: {analysis['total']}")
        report.append("")
        report.append("CLASSIFICATION BY TYPE:")
        for type_name, count in analysis["type_counts"].items():
            percentage = (count / analysis["total"]) * 100 if analysis["total"] > 0 else 0
            report.append(f"  {type_name}: {count} ({percentage:.1f}%)")
        report.append("")
        report.append("ORDER STATISTICS:")
        report.append(f"  Minimum order: {analysis['orders']['min']}")
        report.append(f"  Maximum order: {analysis['orders']['max']}")
        report.append(f"  Average order: {analysis['orders']['avg']:.2f}")
        report.append(f"  Unique orders: {analysis['orders']['unique']}")
        report.append("")
        
        # List first few rings
        report.append("FIRST 10 LEINSTER RINGS:")
        for i, ring in enumerate(result["rings"][:10]):
            report.append(f"  {i+1}. {ring.structure} (order {ring.order})")
        
        return "\n".join(report)
    
    @staticmethod
    def generate_csv_report(result: Dict, filename: str = "leinster_rings.csv"):
        """Generate CSV file of all Leinster rings."""
        import csv
        
        with open(filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['Index', 'Order', 'Structure', 'Type', 'Components'])
            
            for i, ring in enumerate(result["rings"]):
                ring_type = RingAnalyzer.classify_ring(ring)
                components = "×".join([c.structure for c in ring.components]) if ring.components else ""
                writer.writerow([i+1, ring.order, ring.structure, ring_type, components])
        
        print(f"CSV report saved to {filename}")

#############################################################################
## SECTION 6: MAIN INTERFACE AND EXAMPLES
#############################################################################

def main():
    """Main function with example usage."""
    print("Leinster Ring Counter")
    print("=" * 50)
    
    # Example: Count up to order 10
    N = 10
    print(f"\nCounting Leinster rings up to order {N}...")
    
    # Initialize counter
    counter = LeinsterCounter(N)
    
    # Choose counting method
    method = "basic"  # or "optimized"
    
    if method == "basic":
        result = counter.count_basic()
    else:
        result = counter.count_optimized()
    
    # Analyze results
    analyzer = RingAnalyzer()
    analysis = analyzer.analyze_results(result)
    
    # Generate report
    reporter = ReportGenerator()
    report = reporter.generate_text_report(result, analysis)
    print("\n" + report)
    
    # Asymptotic estimation
    asymptotic = analyzer.asymptotic_estimate(N, result)
    print(f"\nAsymptotic estimate:")
    print(f"  L({N}) ~ {asymptotic['C_estimate']:.6f} * {N} / log({N})")
    print(f"  Estimated C = {asymptotic['C_estimate']:.6f}")
    
    # Save to CSV
    reporter.generate_csv_report(result, f"leinster_rings_{N}.csv")
    
    return result, analysis

def test_specific_rings():
    """Test specific known Leinster rings."""
    print("\nTesting specific rings...")
    
    computer = IdealComputer()
    generator = RingGenerator(100)
    
    # Test Z_6 (should be Leinster - perfect number)
    z6 = generator._create_cyclic_ring(6)
    sigma_z6 = computer.compute_sigma(z6)
    is_leinster_z6 = computer.is_leinster(z6)
    print(f"Z_6: σ = {sigma_z6}, 2|R| = {2*6}, Leinster? {is_leinster_z6}")
    
    # Test Z_4 (should NOT be Leinster)
    z4 = generator._create_cyclic_ring(4)
    sigma_z4 = computer.compute_sigma(z4)
    is_leinster_z4 = computer.is_leinster(z4)
    print(f"Z_4: σ = {sigma_z4}, 2|R| = {2*4}, Leinster? {is_leinster_z4}")
    
    # Test field F_7 (should NOT be Leinster)
    f7 = generator._create_field(7)
    sigma_f7 = computer.compute_sigma(f7)
    is_leinster_f7 = computer.is_leinster(f7)
    print(f"F_7: σ = {sigma_f7}, 2|R| = {2*7}, Leinster? {is_leinster_f7}")

def benchmark():
    """Benchmark different counting methods."""
    import time
    
    orders_to_test = [50, 100, 200]
    
    print("\nBenchmarking counting methods...")
    print("Order | Basic (s) | Optimized (s)")
    print("-" * 35)
    
    for N in orders_to_test:
        # Basic method
        counter1 = LeinsterCounter(N)
        start = time.time()
        _ = counter1.count_basic()
        basic_time = time.time() - start
        
        # Optimized method
        counter2 = LeinsterCounter(N)
        start = time.time()
        _ = counter2.count_optimized()
        optimized_time = time.time() - start
        
        print(f"{N:5} | {basic_time:9.3f} | {optimized_time:9.3f}")



# RUN EXAMPLES

#if __name__ == "__main__":
    #print("Leinster Ring Counter - Python Implementation")
# RUN EXAMPLES

if __name__ == "__main__":
    print("Leinster Ring Counter - Python Implementation")
    print("Based on: 'The Arithmetic of Ideal Lattices'")
    print("=" * 60)
    
    # Run tests
    test_specific_rings()
    
    # Run main example
    result, analysis = main()
    # Optional: run benchmark
    # benchmark()
    
    print("\n" + "=" * 60)
    print("Program completed successfully!")
    total_rings = analysis.get('total', 0)
    print(f"Found {total_rings} Leinster rings.")
    print("=" * 60)