# HFSG-zkproof: Zero-Knowledge Proof of Handlebody Spines via Trapping Homomorphisms

This repository implements a **zero-knowledge proof (zk-SNARK)** that a hidden graph, embedded in a public cubical handlebody, is an **essential spine** of the handlebody.  
It combines two lines of work:

1. **Topological foundations** – The paper *Homological Foundations of Spatial Graph Spines and Trapping Homomorphisms* provides a complete graph‑theoretic characterization of *well‑cut* essential spines: a graph can be a spine if and only if every $2$-edge‑connected component is a simple cycle (i.e., a “tree of cycles”).  
2. **ZK circuit & witness** – The paper *A Zero-Knowledge Proof of Handlebody Spine Embedding via Integer Homology* translates the topological criteria (spanning arborescence, integer homology isomorphism, unimodular intersection matrix) into an R1CS circuit written in **Circom 2.1**, accompanied by a **Rust** witness generator that computes all required inputs.

### What does it prove?
The prover convinces the verifier that:
- The hidden graph $G$ lies inside the public handlebody $N$.
- $G$ is connected and has the correct genus $g = \beta_1(G)$.
- The inclusion $G \hookrightarrow N$ induces an isomorphism $H_1(G;\mathbb{Z}) \to H_1(N;\mathbb{Z})$, verified through an integer intersection matrix with determinant $\pm 1$.

By Nielsen’s theorem, this guarantees $G$ is a spine of $N$, without revealing the graph itself.

