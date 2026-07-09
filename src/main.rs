use ark_bn254::Fr;
use ndarray::{Array2, Array1};
use std::collections::VecDeque;

pub struct SpineWitnessGenerator {
    pub v_a: Vec<bool>,
    pub g_e: Vec<bool>,
    pub v_count: usize,
    pub e_count: usize,
    pub genus: usize,
    pub hadamard_bound: i128,
}

impl SpineWitnessGenerator {
    pub fn build_spanning_arborescence(&self, root: usize, edges: &[(usize, usize)]) 
        -> (Vec<[bool; 2]>, Vec<usize>) 
    {
        let mut t_dir = vec![[false; 2]; self.e_count];
        let mut heights = vec![0; self.v_count];
        let mut visited = vec![false; self.v_count];
        
        let mut queue = VecDeque::new();
        queue.push_back((root, 0));
        visited[root] = true;

        while let Some((curr_v, depth)) = queue.pop_front() {
            heights[curr_v] = depth;

            for (edge_idx, &(u, v)) in edges.iter().enumerate() {
                if !self.g_e[edge_idx] { continue; }

                if u == curr_v && !visited[v] {
                    t_dir[edge_idx][0] = true;
                    visited[v] = true;
                    queue.push_back((v, depth + 1));
                } else if v == curr_v && !visited[u] {
                    t_dir[edge_idx][1] = true;
                    visited[u] = true;
                    queue.push_back((u, depth + 1));
                }
            }
        }
        (t_dir, heights)
    }

    pub fn extract_cycle_flows(&self, t_dir: &Vec<[bool; 2]>, edges: &[(usize, usize)]) -> Array2<i32> {
        let mut flows = Array2::zeros((self.e_count, self.genus));
        let mut cycle_idx = 0;

        for (e_idx, &(u, v)) in edges.iter().enumerate() {
            if self.g_e[e_idx] && !t_dir[e_idx][0] && !t_dir[e_idx][1] {
                flows[[e_idx, cycle_idx]] = 1;
                cycle_idx += 1;
                if cycle_idx == self.genus { break; }
            }
        }
        flows
    }

    pub fn compute_inverse_shifted(&self, m_matrix: &Array2<i32>) -> Array2<Fr> {
        let w_matrix_int = integer_matrix_inverse(m_matrix); 

        let mut w_hat = Array2::zeros((self.genus, self.genus));
        for i in 0..self.genus {
            for j in 0..self.genus {
                let shifted_val: i128 = w_matrix_int[[i, j]] as i128 + self.hadamard_bound;
                assert!(shifted_val >= 0);
                assert!(shifted_val <= 2 * self.hadamard_bound);
                w_hat[[i, j]] = Fr::from(shifted_val as u64);
            }
        }
        w_hat
    }
}

fn integer_matrix_inverse(m: &Array2<i32>) -> Array2<i32> {
    let n = m.nrows();
    assert_eq!(n, m.ncols(), "ماتریس باید مربعی باشد");
    if n == 0 {
        return m.clone();
    }
    
    let det = determinant_i32(m);
    assert!(det == 1 || det == -1, "ماتریس unimodular نیست، دترمینان = {}", det);

    let mut adj = Array2::zeros((n, n));
    for i in 0..n {
        for j in 0..n {
            let minor = minor_matrix(m, i, j);
            let cofactor = if (i + j) % 2 == 0 { 1 } else { -1 } * determinant_i32(&minor);
            adj[[j, i]] = cofactor;
        }
    }
    adj.mapv(|x| x * det)
}

fn determinant_i32(m: &Array2<i32>) -> i32 {
    let n = m.nrows();
    if n == 1 {
        return m[[0, 0]];
    }
    if n == 2 {
        return m[[0, 0]] * m[[1, 1]] - m[[0, 1]] * m[[1, 0]];
    }
    let mut det = 0;
    for j in 0..n {
        let minor = minor_matrix(m, 0, j);
        let sign = if j % 2 == 0 { 1 } else { -1 };
        det += sign * m[[0, j]] * determinant_i32(&minor);
    }
    det
}

fn minor_matrix(m: &Array2<i32>, row: usize, col: usize) -> Array2<i32> {
    let n = m.nrows();
    let mut minor = Array2::zeros((n - 1, n - 1));
    let mut r = 0;
    for i in 0..n {
        if i == row { continue; }
        let mut c = 0;
        for j in 0..n {
            if j == col { continue; }
            minor[[r, c]] = m[[i, j]];
            c += 1;
        }
        r += 1;
    }
    minor
}

fn main() {
    let v_a = vec![true; 5];
    let g_e = vec![true; 4];
    let edges = vec![(0,1), (1,2), (2,3), (3,4)];
    let genus = 1;
    let bound = 10;

    let spine_gen = SpineWitnessGenerator {
        v_a,
        g_e: g_e.clone(),
        v_count: 5,
        e_count: 4,
        genus,
        hadamard_bound: bound,
    };

    let (t_dir, heights) = spine_gen.build_spanning_arborescence(0, &edges);
    println!("heights: {:?}", heights);
    println!("t_dir: {:?}", t_dir);

    let m = Array2::from_shape_vec((genus, genus), vec![1]).unwrap();
    let w_hat = spine_gen.compute_inverse_shifted(&m);
    println!("W_hat: {:?}", w_hat);
}