pragma circom 2.1.0;

template Num2Bits(n) {
    signal input in;
    signal output out[n];
    var lc1 = 0;
    var e2 = 1;
    for (var i = 0; i < n; i++) {
        out[i] <-- (in >> i) & 1;
        out[i] * (out[i] - 1) === 0;
        lc1 += out[i] * e2;
        e2 = e2 + e2;
    }
    lc1 === in;
}

template SpanningArborescenceEnforcer(V, E, K, G_GENUS) {
    signal input v_A[V];
    signal input g_e[E];
    signal input t_dir[E][2];
    signal input h[V];
    signal input b_v_k[V][K];
    signal input is_root[V];
    signal input ends[E][2];

    var sum_vA = 0;
    var sum_ge = 0;
    var sum_roots = 0;

    for (var v = 0; v < V; v++) {
        var height_sum = 0;
        var exp = 1;
        for (var k = 0; k < K; k++) {
            b_v_k[v][k] * (1 - b_v_k[v][k]) === 0;
            height_sum += b_v_k[v][k] * exp;
            exp *= 2;
        }
        h[v] === height_sum;

        is_root[v] * (1 - is_root[v]) === 0;
        is_root[v] * h[v] === 0;
        sum_roots += is_root[v];
        sum_vA += v_A[v];
    }
    sum_roots === 1;

    signal incoming[V];
    var in_counts[V];
    for (var v = 0; v < V; v++) in_counts[v] = 0;

    for (var e = 0; e < E; e++) {
        t_dir[e][0] * (1 - t_dir[e][0]) === 0;
        t_dir[e][1] * (1 - t_dir[e][1]) === 0;
        (t_dir[e][0] + t_dir[e][1]) * (1 - g_e[e]) === 0;
        t_dir[e][0] * t_dir[e][1] === 0;
        sum_ge += g_e[e];
    }

    sum_ge - sum_vA + 1 === G_GENUS;
}

template HomologyIsomorphismVerifier(V, E, G_GENUS, B_W, W_BITS) {
    signal input g_e[E];
    signal input x_e_i[E][G_GENUS];
    signal input s_v_e[V][E];
    signal input lambda_e_j[E][G_GENUS];
    signal input W_hat[G_GENUS][G_GENUS];

    signal M[G_GENUS][G_GENUS];

    component x_sq[E][G_GENUS];
    for (var e = 0; e < E; e++) {
        for (var i = 0; i < G_GENUS; i++) {
            x_sq[e][i] = x_e_i[e][i] * x_e_i[e][i];
            x_sq[e][i] * (1 - x_sq[e][i]) === 0;
            x_e_i[e][i] * (1 - g_e[e]) === 0;
        }
    }

    for (var v = 0; v < V; v++) {
        for (var i = 0; i < G_GENUS; i++) {
            var flow_sum = 0;
            for (var e = 0; e < E; e++) {
                flow_sum += s_v_e[v][e] * x_e_i[e][i];
            }
            flow_sum === 0;
        }
    }

    for (var i = 0; i < G_GENUS; i++) {
        for (var j = 0; j < G_GENUS; j++) {
            var m_sum = 0;
            for (var e = 0; e < E; e++) {
                m_sum += x_e_i[e][i] * lambda_e_j[e][j];
            }
            M[i][j] <== m_sum;
        }
    }

    component range_checks[G_GENUS][G_GENUS];
    for (var i = 0; i < G_GENUS; i++) {
        for (var j = 0; j < G_GENUS; j++) {
            range_checks[i][j] = Num2Bits(W_BITS);
            range_checks[i][j].in <== W_hat[i][j];
        }
    }

    for (var i = 0; i < G_GENUS; i++) {
        for (var j = 0; j < G_GENUS; j++) {
            var dot_product = 0;
            for (var k = 0; k < G_GENUS; k++) {
                dot_product += M[i][k] * (W_hat[k][j] - B_W);
            }
            if (i == j) {
                dot_product === 1;
            } else {
                dot_product === 0;
            }
        }
    }
}

template HandlebodySpineProof(V, E, K, G_GENUS, B_W, W_BITS) {
    signal input inN_v[V];
    signal input inN_e[E];
    signal input s_v_e[V][E];
    signal input lambda_e_j[E][G_GENUS];
    signal input ends[E][2];

    signal input v_A[V];
    signal input g_e[E];
    signal input t_dir[E][2];
    signal input h[V];
    signal input b_v_k[V][K];
    signal input is_root[V];
    signal input x_e_i[E][G_GENUS];
    signal input W_hat[G_GENUS][G_GENUS];

    for(var v = 0; v < V; v++){
        v_A[v] * (1 - v_A[v]) === 0;
        v_A[v] * (1 - inN_v[v]) === 0;
    }
    for(var e = 0; e < E; e++){
        g_e[e] * (1 - g_e[e]) === 0;
        g_e[e] * (1 - inN_e[e]) === 0;
    }

    component arb = SpanningArborescenceEnforcer(V, E, K, G_GENUS);
    arb.v_A <== v_A;
    arb.g_e <== g_e;
    arb.t_dir <== t_dir;
    arb.h <== h;
    arb.b_v_k <== b_v_k;
    arb.is_root <== is_root;
    arb.ends <== ends;

    component hom = HomologyIsomorphismVerifier(V, E, G_GENUS, B_W, W_BITS);
    hom.g_e <== g_e;
    hom.x_e_i <== x_e_i;
    hom.s_v_e <== s_v_e;
    hom.lambda_e_j <== lambda_e_j;
    hom.W_hat <== W_hat;
}
