pub type Point = Vec<f64>;

/// K-means clustering (Lloyd’s algorithm).
/// Pseudocode reference: Wikipedia — “K-means clustering”
/// <https://en.wikipedia.org/wiki/K-means_clustering>

pub fn kmeans(points: &[Point], mut centroids: Vec<Point>) -> Vec<Point> {
    if points.is_empty() || centroids.is_empty() {
        return Vec::new();
    }

    let k = centroids.len();
    let mut converged = false;
    let mut iter = 0usize;

    while !converged {
        let mut clusters: Vec<Vec<usize>> = vec![Vec::new(); k];

        for (i, point) in points.iter().enumerate() {
            let mut closest_index = 0usize;
            let mut min_distance = distance(point, &centroids[0]);

            for j in 1..k {
                let d = distance(point, &centroids[j]);
                if d < min_distance {
                    min_distance = d;
                    closest_index = j;
                }
            }
            clusters[closest_index].push(i);
        }

        let mut new_centroids: Vec<Point> = Vec::with_capacity(k);
        for i in 0..k {
            if clusters[i].is_empty() {
                new_centroids.push(centroids[i].clone());
            } else {
                let mut cluster_points: Vec<&Point> = Vec::with_capacity(clusters[i].len());
                for &ix in &clusters[i] {
                    cluster_points.push(&points[ix]);
                }
                new_centroids.push(calculate_centroid(&cluster_points));
            }
        }

        converged = new_centroids == centroids;
        centroids = new_centroids;

        iter += 1;
        if iter > 300 {
            break;
        }
    }

    centroids
}

pub fn calculate_centroid(points_in_cluster: &[&Point]) -> Point {
    let d = points_in_cluster[0].len();
    let mut centroid = vec![0.0; d];

    for p in points_in_cluster {
        for i in 0..d {
            centroid[i] += p[i];
        }
    }
    let n = points_in_cluster.len() as f64;
    for i in 0..d {
        centroid[i] /= n;
    }
    centroid
}

pub fn distance(a: &Point, b: &Point) -> f64 {
    let mut sum = 0.0;
    for i in 0..a.len() {
        let d = a[i] - b[i];
        sum += d * d;
    }
    sum.sqrt()
}
