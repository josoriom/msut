pub type Point = Vec<f64>;

fn dist(a: &Point, b: &Point) -> f64 {
    let mut s = 0.0;
    for i in 0..a.len() {
        let d = a[i] - b[i];
        s += d * d;
    }
    s.sqrt()
}

fn mean(ps: &[&Point]) -> Point {
    let d = ps[0].len();
    let mut m = vec![0.0; d];
    for p in ps {
        for i in 0..d {
            m[i] += p[i];
        }
    }
    let n = ps.len() as f64;
    for i in 0..d {
        m[i] /= n;
    }
    m
}

pub fn kmeans(points: &[Point], mut centroids: Vec<Point>) -> Vec<Point> {
    if points.is_empty() || centroids.is_empty() {
        return Vec::new();
    }
    let k = centroids.len();
    let mut converged = false;
    let mut it = 0usize;

    while !converged {
        let mut groups: Vec<Vec<usize>> = vec![Vec::new(); k];

        for i in 0..points.len() {
            let p = &points[i];
            let mut idx = 0usize;
            let mut best = dist(p, &centroids[0]);
            for j in 1..k {
                let d = dist(p, &centroids[j]);
                if d < best {
                    best = d;
                    idx = j;
                }
            }
            groups[idx].push(i);
        }

        let mut next = Vec::with_capacity(k);
        for gi in 0..k {
            if groups[gi].is_empty() {
                next.push(centroids[gi].clone());
            } else {
                let mut refs: Vec<&Point> = Vec::with_capacity(groups[gi].len());
                for &ix in &groups[gi] {
                    refs.push(&points[ix]);
                }
                next.push(mean(&refs));
            }
        }

        converged = next == centroids;
        centroids = next;
        it += 1;
        if it > 300 {
            break;
        }
    }
    centroids
}
