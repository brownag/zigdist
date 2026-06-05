const std = @import("std");

// Import R headers
const r = @cImport({
    @cInclude("R.h");
    @cInclude("Rinternals.h");
});

fn clamp_threads(requested: i32) usize {
    const system_cores = std.Thread.getCpuCount() catch 1;
    var count: usize = if (requested > 0) @as(usize, @intCast(requested)) else 1;
    if (count > system_cores) count = system_cores;
    return count;
}

fn transpose_matrix(src: [*]const f64, dst: []f64, nrows: usize, ncols: usize) void {
    var i: usize = 0;
    while (i < nrows) : (i += 1) {
        var j: usize = 0;
        while (j < ncols) : (j += 1) {
            dst[i * ncols + j] = src[j * nrows + i];
        }
    }
}

fn transpose_int_matrix(src: [*]const i32, dst: []i32, nrows: usize, ncols: usize) void {
    var i: usize = 0;
    while (i < nrows) : (i += 1) {
        var j: usize = 0;
        while (j < ncols) : (j += 1) {
            dst[i * ncols + j] = src[j * nrows + i];
        }
    }
}

// Euclidean Worker Functions (no NA)

fn euclidean_self_worker_fast(
    n: usize,
    p: usize,
    transposed: [*]const f64,
    res: [*]f64,
    thread_idx: usize,
    thread_count: usize,
) void {
    var idx1: usize = thread_idx;
    while (idx1 < n) : (idx1 += thread_count) {
        res[idx1 * n + idx1] = 0.0;
        var idx2: usize = idx1 + 1;
        while (idx2 < n) : (idx2 += 1) {
            var sum: f64 = 0.0;
            var j: usize = 0;
            const r1 = idx1 * p;
            const r2 = idx2 * p;
            while (j < p) : (j += 1) {
                const diff = transposed[r1 + j] - transposed[r2 + j];
                sum += diff * diff;
            }
            const dist = @sqrt(sum);
            res[idx2 * n + idx1] = dist;
            res[idx1 * n + idx2] = dist;
        }
    }
}

fn euclidean_cross_worker_fast(
    nx: usize,
    ny: usize,
    p: usize,
    x_trans: [*]const f64,
    y_trans: [*]const f64,
    res: [*]f64,
    thread_idx: usize,
    thread_count: usize,
) void {
    var idx2: usize = thread_idx;
    while (idx2 < ny) : (idx2 += thread_count) {
        var idx1: usize = 0;
        const r2 = idx2 * p;
        while (idx1 < nx) : (idx1 += 1) {
            var sum: f64 = 0.0;
            var j: usize = 0;
            const r1 = idx1 * p;
            while (j < p) : (j += 1) {
                const diff = x_trans[r1 + j] - y_trans[r2 + j];
                sum += diff * diff;
            }
            res[idx2 * nx + idx1] = @sqrt(sum);
        }
    }
}

// Euclidean Worker Functions (with NA)

fn euclidean_self_worker_na(
    n: usize,
    p: usize,
    transposed: [*]const f64,
    res: [*]f64,
    thread_idx: usize,
    thread_count: usize,
) void {
    var idx1: usize = thread_idx;
    while (idx1 < n) : (idx1 += thread_count) {
        res[idx1 * n + idx1] = 0.0;
        var idx2: usize = idx1 + 1;
        while (idx2 < n) : (idx2 += 1) {
            var sum: f64 = 0.0;
            var present_count: f64 = 0.0;
            var j: usize = 0;
            const r1 = idx1 * p;
            const r2 = idx2 * p;
            while (j < p) : (j += 1) {
                const v1 = transposed[r1 + j];
                const v2 = transposed[r2 + j];
                // Check for NaN (missing)
                if (!std.math.isNan(v1) and !std.math.isNan(v2)) {
                    const diff = v1 - v2;
                    sum += diff * diff;
                    present_count += 1.0;
                }
            }

            if (present_count > 0.0) {
                const scale = @as(f64, @floatFromInt(p)) / present_count;
                const dist = @sqrt(sum * scale);
                res[idx2 * n + idx1] = dist;
                res[idx1 * n + idx2] = dist;
            } else {
                res[idx2 * n + idx1] = std.math.nan(f64);
                res[idx1 * n + idx2] = std.math.nan(f64);
            }
        }
    }
}

fn euclidean_cross_worker_na(
    nx: usize,
    ny: usize,
    p: usize,
    x_trans: [*]const f64,
    y_trans: [*]const f64,
    res: [*]f64,
    thread_idx: usize,
    thread_count: usize,
) void {
    var idx2: usize = thread_idx;
    while (idx2 < ny) : (idx2 += thread_count) {
        var idx1: usize = 0;
        const r2 = idx2 * p;
        while (idx1 < nx) : (idx1 += 1) {
            var sum: f64 = 0.0;
            var present_count: f64 = 0.0;
            var j: usize = 0;
            const r1 = idx1 * p;
            while (j < p) : (j += 1) {
                const v1 = x_trans[r1 + j];
                const v2 = y_trans[r2 + j];
                if (!std.math.isNan(v1) and !std.math.isNan(v2)) {
                    const diff = v1 - v2;
                    sum += diff * diff;
                    present_count += 1.0;
                }
            }

            if (present_count > 0.0) {
                const scale = @as(f64, @floatFromInt(p)) / present_count;
                res[idx2 * nx + idx1] = @sqrt(sum * scale);
            } else {
                res[idx2 * nx + idx1] = std.math.nan(f64);
            }
        }
    }
}

// Gower Worker Functions (no NA)

fn gower_self_worker_fast(
    n: usize,
    p_num: usize,
    p_cat: usize,
    x_num_trans: [*]const f64,
    x_cat_trans: ?[*]const i32,
    ranges: [*]const f64,
    w_num: [*]const f64,
    w_cat: ?[*]const f64,
    res: [*]f64,
    thread_idx: usize,
    thread_count: usize,
) void {
    var idx1: usize = thread_idx;
    while (idx1 < n) : (idx1 += thread_count) {
        res[idx1 * n + idx1] = 0.0;
        var idx2: usize = idx1 + 1;
        while (idx2 < n) : (idx2 += 1) {
            var weighted_sum: f64 = 0.0;
            var total_weight: f64 = 0.0;

            // Numeric features
            var k: usize = 0;
            const r1_num = idx1 * p_num;
            const r2_num = idx2 * p_num;
            while (k < p_num) : (k += 1) {
                const rng = ranges[k];
                if (rng > 0.0) {
                    const diff = x_num_trans[r1_num + k] - x_num_trans[r2_num + k];
                    const abs_diff = if (diff < 0.0) -diff else diff;
                    weighted_sum += w_num[k] * (abs_diff / rng);
                    total_weight += w_num[k];
                }
            }

            // Categorical features
            if (p_cat > 0) {
                const cat_ptr = x_cat_trans.?;
                const w_cat_ptr = w_cat.?;
                k = 0;
                const r1_cat = idx1 * p_cat;
                const r2_cat = idx2 * p_cat;
                while (k < p_cat) : (k += 1) {
                    const diff: f64 = if (cat_ptr[r1_cat + k] != cat_ptr[r2_cat + k]) 1.0 else 0.0;
                    weighted_sum += w_cat_ptr[k] * diff;
                    total_weight += w_cat_ptr[k];
                }
            }

            const dist = if (total_weight > 0.0) weighted_sum / total_weight else std.math.nan(f64);
            res[idx2 * n + idx1] = dist;
            res[idx1 * n + idx2] = dist;
        }
    }
}

fn gower_cross_worker_fast(
    nx: usize,
    ny: usize,
    p_num: usize,
    p_cat: usize,
    x_num_trans: [*]const f64,
    x_cat_trans: ?[*]const i32,
    y_num_trans: [*]const f64,
    y_cat_trans: ?[*]const i32,
    ranges: [*]const f64,
    w_num: [*]const f64,
    w_cat: ?[*]const f64,
    res: [*]f64,
    thread_idx: usize,
    thread_count: usize,
) void {
    var idx2: usize = thread_idx;
    while (idx2 < ny) : (idx2 += thread_count) {
        var idx1: usize = 0;
        const r2_num = idx2 * p_num;
        const r2_cat = idx2 * p_cat;
        while (idx1 < nx) : (idx1 += 1) {
            var weighted_sum: f64 = 0.0;
            var total_weight: f64 = 0.0;

            var k: usize = 0;
            const r1_num = idx1 * p_num;
            while (k < p_num) : (k += 1) {
                const rng = ranges[k];
                if (rng > 0.0) {
                    const diff = x_num_trans[r1_num + k] - y_num_trans[r2_num + k];
                    const abs_diff = if (diff < 0.0) -diff else diff;
                    weighted_sum += w_num[k] * (abs_diff / rng);
                    total_weight += w_num[k];
                }
            }

            if (p_cat > 0) {
                const x_cat_ptr = x_cat_trans.?;
                const y_cat_ptr = y_cat_trans.?;
                const w_cat_ptr = w_cat.?;
                k = 0;
                const r1_cat = idx1 * p_cat;
                while (k < p_cat) : (k += 1) {
                    const diff: f64 = if (x_cat_ptr[r1_cat + k] != y_cat_ptr[r2_cat + k]) 1.0 else 0.0;
                    weighted_sum += w_cat_ptr[k] * diff;
                    total_weight += w_cat_ptr[k];
                }
            }

            res[idx2 * nx + idx1] = if (total_weight > 0.0) weighted_sum / total_weight else std.math.nan(f64);
        }
    }
}

// Gower Worker Functions (with NA)

fn gower_self_worker_na(
    n: usize,
    p_num: usize,
    p_cat: usize,
    x_num_trans: [*]const f64,
    x_cat_trans: ?[*]const i32,
    ranges: [*]const f64,
    w_num: [*]const f64,
    w_cat: ?[*]const f64,
    res: [*]f64,
    thread_idx: usize,
    thread_count: usize,
) void {
    var idx1: usize = thread_idx;
    while (idx1 < n) : (idx1 += thread_count) {
        res[idx1 * n + idx1] = 0.0;
        var idx2: usize = idx1 + 1;
        while (idx2 < n) : (idx2 += 1) {
            var weighted_sum: f64 = 0.0;
            var total_weight: f64 = 0.0;

            // Numeric features
            var k: usize = 0;
            const r1_num = idx1 * p_num;
            const r2_num = idx2 * p_num;
            while (k < p_num) : (k += 1) {
                const rng = ranges[k];
                if (rng > 0.0) {
                    const v1 = x_num_trans[r1_num + k];
                    const v2 = x_num_trans[r2_num + k];
                    if (!std.math.isNan(v1) and !std.math.isNan(v2)) {
                        const diff = v1 - v2;
                        const abs_diff = if (diff < 0.0) -diff else diff;
                        weighted_sum += w_num[k] * (abs_diff / rng);
                        total_weight += w_num[k];
                    }
                }
            }

            // Categorical features
            if (p_cat > 0) {
                const cat_ptr = x_cat_trans.?;
                const w_cat_ptr = w_cat.?;
                k = 0;
                const r1_cat = idx1 * p_cat;
                const r2_cat = idx2 * p_cat;
                while (k < p_cat) : (k += 1) {
                    const val1 = cat_ptr[r1_cat + k];
                    const val2 = cat_ptr[r2_cat + k];
                    // Check for NA_INTEGER (-2147483648)
                    if (val1 != -2147483648 and val2 != -2147483648) {
                        const diff: f64 = if (val1 != val2) 1.0 else 0.0;
                        weighted_sum += w_cat_ptr[k] * diff;
                        total_weight += w_cat_ptr[k];
                    }
                }
            }

            const dist = if (total_weight > 0.0) weighted_sum / total_weight else std.math.nan(f64);
            res[idx2 * n + idx1] = dist;
            res[idx1 * n + idx2] = dist;
        }
    }
}

fn gower_cross_worker_na(
    nx: usize,
    ny: usize,
    p_num: usize,
    p_cat: usize,
    x_num_trans: [*]const f64,
    x_cat_trans: ?[*]const i32,
    y_num_trans: [*]const f64,
    y_cat_trans: ?[*]const i32,
    ranges: [*]const f64,
    w_num: [*]const f64,
    w_cat: ?[*]const f64,
    res: [*]f64,
    thread_idx: usize,
    thread_count: usize,
) void {
    var idx2: usize = thread_idx;
    while (idx2 < ny) : (idx2 += thread_count) {
        var idx1: usize = 0;
        const r2_num = idx2 * p_num;
        const r2_cat = idx2 * p_cat;
        while (idx1 < nx) : (idx1 += 1) {
            var weighted_sum: f64 = 0.0;
            var total_weight: f64 = 0.0;

            // Numeric features
            var k: usize = 0;
            const r1_num = idx1 * p_num;
            while (k < p_num) : (k += 1) {
                const rng = ranges[k];
                if (rng > 0.0) {
                    const v1 = x_num_trans[r1_num + k];
                    const v2 = y_num_trans[r2_num + k];
                    if (!std.math.isNan(v1) and !std.math.isNan(v2)) {
                        const diff = v1 - v2;
                        const abs_diff = if (diff < 0.0) -diff else diff;
                        weighted_sum += w_num[k] * (abs_diff / rng);
                        total_weight += w_num[k];
                    }
                }
            }

            // Categorical features
            if (p_cat > 0) {
                const x_cat_ptr = x_cat_trans.?;
                const y_cat_ptr = y_cat_trans.?;
                const w_cat_ptr = w_cat.?;
                k = 0;
                const r1_cat = idx1 * p_cat;
                while (k < p_cat) : (k += 1) {
                    const val1 = x_cat_ptr[r1_cat + k];
                    const val2 = y_cat_ptr[r2_cat + k];
                    if (val1 != -2147483648 and val2 != -2147483648) {
                        const diff: f64 = if (val1 != val2) 1.0 else 0.0;
                        weighted_sum += w_cat_ptr[k] * diff;
                        total_weight += w_cat_ptr[k];
                    }
                }
            }

            res[idx2 * nx + idx1] = if (total_weight > 0.0) weighted_sum / total_weight else std.math.nan(f64);
        }
    }
}

// --- Thread spawning helpers ---

fn spawn_euclidean_self_workers(
    thread_count: usize,
    small_threshold: usize,
    has_na: bool,
    worker_na: fn(usize, usize, [*]const f64, [*]f64, usize, usize) void,
    worker_fast: fn(usize, usize, [*]const f64, [*]f64, usize, usize) void,
    n: usize,
    p: usize,
    transposed: [*]const f64,
    res: [*]f64,
) void {
    if (thread_count <= 1 or n < small_threshold) {
        if (has_na) {
            worker_na(n, p, transposed, res, 0, 1);
        } else {
            worker_fast(n, p, transposed, res, 0, 1);
        }
        return;
    }

    var threads = std.heap.page_allocator.alloc(std.Thread, thread_count) catch {
        if (has_na) {
            worker_na(n, p, transposed, res, 0, 1);
        } else {
            worker_fast(n, p, transposed, res, 0, 1);
        }
        return;
    };
    defer std.heap.page_allocator.free(threads);

    var t: usize = 0;
    while (t < thread_count) : (t += 1) {
        if (has_na) {
            threads[t] = std.Thread.spawn(.{}, worker_na, .{ n, p, transposed, res, t, thread_count }) catch {
                worker_na(n, p, transposed, res, t, thread_count);
                threads = threads[0..t];
                break;
            };
        } else {
            threads[t] = std.Thread.spawn(.{}, worker_fast, .{ n, p, transposed, res, t, thread_count }) catch {
                worker_fast(n, p, transposed, res, t, thread_count);
                threads = threads[0..t];
                break;
            };
        }
    }

    for (threads) |th| th.join();
}

fn spawn_euclidean_cross_workers(
    thread_count: usize,
    small_threshold: usize,
    has_na: bool,
    worker_na: fn(usize, usize, usize, [*]const f64, [*]const f64, [*]f64, usize, usize) void,
    worker_fast: fn(usize, usize, usize, [*]const f64, [*]const f64, [*]f64, usize, usize) void,
    nx: usize,
    ny: usize,
    p: usize,
    x_trans: [*]const f64,
    y_trans: [*]const f64,
    res: [*]f64,
) void {
    if (thread_count <= 1 or ny < small_threshold) {
        if (has_na) {
            worker_na(nx, ny, p, x_trans, y_trans, res, 0, 1);
        } else {
            worker_fast(nx, ny, p, x_trans, y_trans, res, 0, 1);
        }
        return;
    }

    var threads = std.heap.page_allocator.alloc(std.Thread, thread_count) catch {
        if (has_na) {
            worker_na(nx, ny, p, x_trans, y_trans, res, 0, 1);
        } else {
            worker_fast(nx, ny, p, x_trans, y_trans, res, 0, 1);
        }
        return;
    };
    defer std.heap.page_allocator.free(threads);

    var t: usize = 0;
    while (t < thread_count) : (t += 1) {
        if (has_na) {
            threads[t] = std.Thread.spawn(.{}, worker_na, .{ nx, ny, p, x_trans, y_trans, res, t, thread_count }) catch {
                worker_na(nx, ny, p, x_trans, y_trans, res, t, thread_count);
                threads = threads[0..t];
                break;
            };
        } else {
            threads[t] = std.Thread.spawn(.{}, worker_fast, .{ nx, ny, p, x_trans, y_trans, res, t, thread_count }) catch {
                worker_fast(nx, ny, p, x_trans, y_trans, res, t, thread_count);
                threads = threads[0..t];
                break;
            };
        }
    }

    for (threads) |th| th.join();
}

fn spawn_gower_self_workers(
    thread_count: usize,
    small_threshold: usize,
    has_na: bool,
    worker_na: fn(usize, usize, usize, [*]const f64, ?[*]const i32, [*]const f64, [*]const f64, ?[*]const f64, [*]f64, usize, usize) void,
    worker_fast: fn(usize, usize, usize, [*]const f64, ?[*]const i32, [*]const f64, [*]const f64, ?[*]const f64, [*]f64, usize, usize) void,
    n: usize,
    p_num: usize,
    p_cat: usize,
    x_num_trans: [*]const f64,
    x_cat_trans: ?[*]const i32,
    ranges: [*]const f64,
    w_num: [*]const f64,
    w_cat: ?[*]const f64,
    res: [*]f64,
) void {
    if (thread_count <= 1 or n < small_threshold) {
        if (has_na) {
            worker_na(n, p_num, p_cat, x_num_trans, x_cat_trans, ranges, w_num, w_cat, res, 0, 1);
        } else {
            worker_fast(n, p_num, p_cat, x_num_trans, x_cat_trans, ranges, w_num, w_cat, res, 0, 1);
        }
        return;
    }

    var threads = std.heap.page_allocator.alloc(std.Thread, thread_count) catch {
        if (has_na) {
            worker_na(n, p_num, p_cat, x_num_trans, x_cat_trans, ranges, w_num, w_cat, res, 0, 1);
        } else {
            worker_fast(n, p_num, p_cat, x_num_trans, x_cat_trans, ranges, w_num, w_cat, res, 0, 1);
        }
        return;
    };
    defer std.heap.page_allocator.free(threads);

    var t: usize = 0;
    while (t < thread_count) : (t += 1) {
        if (has_na) {
            threads[t] = std.Thread.spawn(.{}, worker_na, .{ n, p_num, p_cat, x_num_trans, x_cat_trans, ranges, w_num, w_cat, res, t, thread_count }) catch {
                worker_na(n, p_num, p_cat, x_num_trans, x_cat_trans, ranges, w_num, w_cat, res, t, thread_count);
                threads = threads[0..t];
                break;
            };
        } else {
            threads[t] = std.Thread.spawn(.{}, worker_fast, .{ n, p_num, p_cat, x_num_trans, x_cat_trans, ranges, w_num, w_cat, res, t, thread_count }) catch {
                worker_fast(n, p_num, p_cat, x_num_trans, x_cat_trans, ranges, w_num, w_cat, res, t, thread_count);
                threads = threads[0..t];
                break;
            };
        }
    }

    for (threads) |th| th.join();
}

fn spawn_gower_cross_workers(
    thread_count: usize,
    small_threshold: usize,
    has_na: bool,
    worker_na: fn(usize, usize, usize, usize, [*]const f64, ?[*]const i32, [*]const f64, ?[*]const i32, [*]const f64, [*]const f64, ?[*]const f64, [*]f64, usize, usize) void,
    worker_fast: fn(usize, usize, usize, usize, [*]const f64, ?[*]const i32, [*]const f64, ?[*]const i32, [*]const f64, [*]const f64, ?[*]const f64, [*]f64, usize, usize) void,
    nx: usize,
    ny: usize,
    p_num: usize,
    p_cat: usize,
    x_num_trans: [*]const f64,
    x_cat_trans: ?[*]const i32,
    y_num_trans: [*]const f64,
    y_cat_trans: ?[*]const i32,
    ranges: [*]const f64,
    w_num: [*]const f64,
    w_cat: ?[*]const f64,
    res: [*]f64,
) void {
    if (thread_count <= 1 or ny < small_threshold) {
        if (has_na) {
            worker_na(nx, ny, p_num, p_cat, x_num_trans, x_cat_trans, y_num_trans, y_cat_trans, ranges, w_num, w_cat, res, 0, 1);
        } else {
            worker_fast(nx, ny, p_num, p_cat, x_num_trans, x_cat_trans, y_num_trans, y_cat_trans, ranges, w_num, w_cat, res, 0, 1);
        }
        return;
    }

    var threads = std.heap.page_allocator.alloc(std.Thread, thread_count) catch {
        if (has_na) {
            worker_na(nx, ny, p_num, p_cat, x_num_trans, x_cat_trans, y_num_trans, y_cat_trans, ranges, w_num, w_cat, res, 0, 1);
        } else {
            worker_fast(nx, ny, p_num, p_cat, x_num_trans, x_cat_trans, y_num_trans, y_cat_trans, ranges, w_num, w_cat, res, 0, 1);
        }
        return;
    };
    defer std.heap.page_allocator.free(threads);

    var t: usize = 0;
    while (t < thread_count) : (t += 1) {
        if (has_na) {
            threads[t] = std.Thread.spawn(.{}, worker_na, .{ nx, ny, p_num, p_cat, x_num_trans, x_cat_trans, y_num_trans, y_cat_trans, ranges, w_num, w_cat, res, t, thread_count }) catch {
                worker_na(nx, ny, p_num, p_cat, x_num_trans, x_cat_trans, y_num_trans, y_cat_trans, ranges, w_num, w_cat, res, t, thread_count);
                threads = threads[0..t];
                break;
            };
        } else {
            threads[t] = std.Thread.spawn(.{}, worker_fast, .{ nx, ny, p_num, p_cat, x_num_trans, x_cat_trans, y_num_trans, y_cat_trans, ranges, w_num, w_cat, res, t, thread_count }) catch {
                worker_fast(nx, ny, p_num, p_cat, x_num_trans, x_cat_trans, y_num_trans, y_cat_trans, ranges, w_num, w_cat, res, t, thread_count);
                threads = threads[0..t];
                break;
            };
        }
    }

    for (threads) |th| th.join();
}

// R API Entry Points

export fn zd_euclidean_dist_self_R(x_sexp: r.SEXP, num_threads_sexp: r.SEXP, has_na_sexp: r.SEXP) r.SEXP {
    const N = r.Rf_nrows(x_sexp);
    const P = r.Rf_ncols(x_sexp);
    const x = r.REAL(x_sexp).?;

    const res_sexp = r.Rf_protect(r.Rf_allocMatrix(r.REALSXP, N, N));
    const res = r.REAL(res_sexp).?;

    const n = @as(usize, @intCast(N));
    const p = @as(usize, @intCast(P));

    const thread_count = clamp_threads(r.INTEGER(num_threads_sexp)[0]);
    const has_na = r.LOGICAL(has_na_sexp)[0] != 0;

    // Allocate and transpose to contiguous row-major buffer
    const transposed = std.heap.page_allocator.alloc(f64, n * p) catch {
        r.Rf_error("Memory allocation failed for transposed matrix in zigdist");
        unreachable;
    };
    defer std.heap.page_allocator.free(transposed);

    transpose_matrix(x, transposed, n, p);

    spawn_euclidean_self_workers(thread_count, 20, has_na,
        euclidean_self_worker_na, euclidean_self_worker_fast,
        n, p, transposed.ptr, res);

    r.Rf_unprotect(1);
    return res_sexp;
}

export fn zd_euclidean_dist_cross_R(
    x_sexp: r.SEXP,
    y_sexp: r.SEXP,
    num_threads_sexp: r.SEXP,
    has_na_sexp: r.SEXP,
) r.SEXP {
    const N_x = r.Rf_nrows(x_sexp);
    const N_y = r.Rf_nrows(y_sexp);
    const P = r.Rf_ncols(x_sexp);
    const x = r.REAL(x_sexp).?;
    const y = r.REAL(y_sexp).?;

    const res_sexp = r.Rf_protect(r.Rf_allocMatrix(r.REALSXP, N_x, N_y));
    const res = r.REAL(res_sexp).?;

    const nx = @as(usize, @intCast(N_x));
    const ny = @as(usize, @intCast(N_y));
    const p = @as(usize, @intCast(P));

    const thread_count = clamp_threads(r.INTEGER(num_threads_sexp)[0]);
    const has_na = r.LOGICAL(has_na_sexp)[0] != 0;

    // Allocate and transpose x
    const x_transposed = std.heap.page_allocator.alloc(f64, nx * p) catch {
        r.Rf_error("Memory allocation failed for x_transposed in zigdist");
        unreachable;
    };
    defer std.heap.page_allocator.free(x_transposed);

    transpose_matrix(x, x_transposed, nx, p);

    // Allocate and transpose y
    const y_transposed = std.heap.page_allocator.alloc(f64, ny * p) catch {
        r.Rf_error("Memory allocation failed for y_transposed in zigdist");
        unreachable;
    };
    defer std.heap.page_allocator.free(y_transposed);

    transpose_matrix(y, y_transposed, ny, p);

    spawn_euclidean_cross_workers(thread_count, 10, has_na,
        euclidean_cross_worker_na, euclidean_cross_worker_fast,
        nx, ny, p, x_transposed.ptr, y_transposed.ptr, res);

    r.Rf_unprotect(1);
    return res_sexp;
}

export fn zd_gower_dist_self_R(
    x_num_sexp: r.SEXP,
    x_cat_sexp: r.SEXP,
    ranges_sexp: r.SEXP,
    w_num_sexp: r.SEXP,
    w_cat_sexp: r.SEXP,
    num_threads_sexp: r.SEXP,
    has_na_sexp: r.SEXP,
) r.SEXP {
    const N = r.Rf_nrows(x_num_sexp);
    const P_num = r.Rf_ncols(x_num_sexp);
    const P_cat = r.Rf_ncols(x_cat_sexp);

    const x_num = r.REAL(x_num_sexp).?;
    const x_cat = if (P_cat > 0) r.INTEGER(x_cat_sexp).? else null;
    const ranges = r.REAL(ranges_sexp).?;
    const w_num = r.REAL(w_num_sexp).?;
    const w_cat = if (P_cat > 0) r.REAL(w_cat_sexp).? else null;

    const res_sexp = r.Rf_protect(r.Rf_allocMatrix(r.REALSXP, N, N));
    const res = r.REAL(res_sexp).?;

    const n = @as(usize, @intCast(N));
    const p_num = @as(usize, @intCast(P_num));
    const p_cat = @as(usize, @intCast(P_cat));

    const thread_count = clamp_threads(r.INTEGER(num_threads_sexp)[0]);
    const has_na = r.LOGICAL(has_na_sexp)[0] != 0;

    // Allocate and transpose numeric features
    const x_num_trans = std.heap.page_allocator.alloc(f64, n * p_num) catch {
        r.Rf_error("Memory allocation failed for x_num_trans in zigdist");
        unreachable;
    };
    defer std.heap.page_allocator.free(x_num_trans);

    transpose_matrix(x_num, x_num_trans, n, p_num);

    // Allocate and transpose categorical features if any
    var x_cat_trans: ?[]i32 = null;
    if (p_cat > 0) {
        x_cat_trans = std.heap.page_allocator.alloc(i32, n * p_cat) catch {
            r.Rf_error("Memory allocation failed for x_cat_trans in zigdist");
            unreachable;
        };
    }
    defer if (x_cat_trans) |slice| std.heap.page_allocator.free(slice);

    if (p_cat > 0) {
        const cat_ptr = x_cat.?;
        const trans_ptr = x_cat_trans.?.ptr;
        transpose_int_matrix(cat_ptr, trans_ptr[0..n * p_cat], n, p_cat);
    }

    const x_cat_trans_ptr = if (x_cat_trans) |slice| slice.ptr else null;

    spawn_gower_self_workers(thread_count, 20, has_na,
        gower_self_worker_na, gower_self_worker_fast,
        n, p_num, p_cat, x_num_trans.ptr, x_cat_trans_ptr, ranges, w_num, w_cat, res);

    r.Rf_unprotect(1);
    return res_sexp;
}

export fn zd_gower_dist_cross_R(
    x_num_sexp: r.SEXP,
    x_cat_sexp: r.SEXP,
    y_num_sexp: r.SEXP,
    y_cat_sexp: r.SEXP,
    ranges_sexp: r.SEXP,
    w_num_sexp: r.SEXP,
    w_cat_sexp: r.SEXP,
    num_threads_sexp: r.SEXP,
    has_na_sexp: r.SEXP,
) r.SEXP {
    const N_x = r.Rf_nrows(x_num_sexp);
    const N_y = r.Rf_nrows(y_num_sexp);
    const P_num = r.Rf_ncols(x_num_sexp);
    const P_cat = r.Rf_ncols(x_cat_sexp);

    const x_num = r.REAL(x_num_sexp).?;
    const x_cat = if (P_cat > 0) r.INTEGER(x_cat_sexp).? else null;
    const y_num = r.REAL(y_num_sexp).?;
    const y_cat = if (P_cat > 0) r.INTEGER(y_cat_sexp).? else null;
    const ranges = r.REAL(ranges_sexp).?;
    const w_num = r.REAL(w_num_sexp).?;
    const w_cat = if (P_cat > 0) r.REAL(w_cat_sexp).? else null;

    const res_sexp = r.Rf_protect(r.Rf_allocMatrix(r.REALSXP, N_x, N_y));
    const res = r.REAL(res_sexp).?;

    const nx = @as(usize, @intCast(N_x));
    const ny = @as(usize, @intCast(N_y));
    const p_num = @as(usize, @intCast(P_num));
    const p_cat = @as(usize, @intCast(P_cat));

    const thread_count = clamp_threads(r.INTEGER(num_threads_sexp)[0]);
    const has_na = r.LOGICAL(has_na_sexp)[0] != 0;

    // Transpose x_num
    const x_num_trans = std.heap.page_allocator.alloc(f64, nx * p_num) catch {
        r.Rf_error("Memory allocation failed for x_num_trans in zigdist");
        unreachable;
    };
    defer std.heap.page_allocator.free(x_num_trans);

    transpose_matrix(x_num, x_num_trans, nx, p_num);

    // Transpose y_num
    const y_num_trans = std.heap.page_allocator.alloc(f64, ny * p_num) catch {
        r.Rf_error("Memory allocation failed for y_num_trans in zigdist");
        unreachable;
    };
    defer std.heap.page_allocator.free(y_num_trans);

    transpose_matrix(y_num, y_num_trans, ny, p_num);

    // Transpose x_cat
    var x_cat_trans: ?[]i32 = null;
    if (p_cat > 0) {
        x_cat_trans = std.heap.page_allocator.alloc(i32, nx * p_cat) catch {
            r.Rf_error("Memory allocation failed for x_cat_trans in zigdist");
            unreachable;
        };
    }
    defer if (x_cat_trans) |slice| std.heap.page_allocator.free(slice);

    if (p_cat > 0) {
        const cat_ptr = x_cat.?;
        const trans_ptr = x_cat_trans.?.ptr;
        transpose_int_matrix(cat_ptr, trans_ptr[0..nx * p_cat], nx, p_cat);
    }

    // Transpose y_cat
    var y_cat_trans: ?[]i32 = null;
    if (p_cat > 0) {
        y_cat_trans = std.heap.page_allocator.alloc(i32, ny * p_cat) catch {
            r.Rf_error("Memory allocation failed for y_cat_trans in zigdist");
            unreachable;
        };
    }
    defer if (y_cat_trans) |slice| std.heap.page_allocator.free(slice);

    if (p_cat > 0) {
        const cat_ptr = y_cat.?;
        const trans_ptr = y_cat_trans.?.ptr;
        transpose_int_matrix(cat_ptr, trans_ptr[0..ny * p_cat], ny, p_cat);
    }

    const x_cat_trans_ptr = if (x_cat_trans) |slice| slice.ptr else null;
    const y_cat_trans_ptr = if (y_cat_trans) |slice| slice.ptr else null;

    spawn_gower_cross_workers(thread_count, 10, has_na,
        gower_cross_worker_na, gower_cross_worker_fast,
        nx, ny, p_num, p_cat, x_num_trans.ptr, x_cat_trans_ptr, y_num_trans.ptr, y_cat_trans_ptr, ranges, w_num, w_cat, res);

    r.Rf_unprotect(1);
    return res_sexp;
}
