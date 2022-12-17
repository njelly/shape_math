use bevy::math::Vec2;
use std::f32::consts::PI;

pub fn aabb_contains_point(min: Vec2, max: Vec2, p: Vec2) -> bool {
    p.x >= min.x && p.x <= max.x && p.y >= min.y && p.y <= max.y
}

pub fn aabb_intersects_aabb(min_a: Vec2, max_a: Vec2, min_b: Vec2, max_b: Vec2) -> bool {
    max_a.x - min_a.x + max_b.x - min_b.x >= max_a.x.max(max_b.x) - min_a.x.min(min_b.x)
        && max_a.y - min_a.y + max_b.y - min_b.y >= max_a.y.max(max_b.y) - min_a.y.min(min_b.y)
}

pub fn get_vertices_aabb(min: Vec2, max: Vec2) -> Vec<Vec2> {
    vec![max, Vec2::new(max.x, min.y), min, Vec2::new(min.x, max.y)]
}

pub fn get_bounding_aabb(points: &Vec<Vec2>) -> (Vec2, Vec2) {
    let mut min = Vec2::new(f32::MAX, f32::MAX);
    let mut max = Vec2::new(f32::MIN, f32::MIN);
    for p in points {
        if p.x < min.x {
            min.x = p.x;
        } else if p.x > max.x {
            max.x = p.x;
        }

        if p.y < min.y {
            min.y = p.y;
        } else if p.y > max.y {
            max.y = p.y;
        }
    }

    (min, max)
}

pub fn circle_contains_point(center: Vec2, radius: f32, p: Vec2) -> bool {
    (p - center).length_squared() <= radius * radius
}

pub fn circle_intersects_aabb(center: Vec2, r: f32, aabb_min: Vec2, aabb_max: Vec2) -> bool {
    if aabb_contains_point(aabb_min, aabb_max, center) {
        return true;
    }

    let vertices = get_vertices_aabb(aabb_min, aabb_max);
    for (i, v) in vertices.iter().enumerate() {
        let next_index = (i + 1) % vertices.len();
        let closest_point_on_edge = closest_point_on_line_segment(*v, vertices[next_index], center);
        if circle_contains_point(center, r, closest_point_on_edge) {
            return true;
        }
    }

    false
}

pub fn circle_intersects_circle(center_a: &Vec2, r_a: &f32, center_b: &Vec2, r_b: &f32) -> bool {
    let r_sum = r_a + r_b;
    (*center_b - *center_a).length_squared() <= r_sum * r_sum
}

pub fn closest_point_on_line(a: Vec2, b: Vec2, p: Vec2) -> Vec2 {
    let ap = p - a;
    let ab = b - a;
    let dist = ap.dot(ab) / ab.length_squared();
    a + ab * dist
}

pub fn closest_point_on_line_segment(a: Vec2, b: Vec2, p: Vec2) -> Vec2 {
    let ap = p - a;
    let ab = b - a;
    let dist = ap.dot(ab) / ab.length_squared();
    if dist <= 0. {
        a
    } else if dist >= 1. {
        b
    } else {
        a + ab * dist
    }
}

pub fn point_is_on_left_side_of_line(a: Vec2, b: Vec2, p: Vec2) -> bool {
    (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x) > 0.
}

pub fn get_longest_edge_of_polygon(vertices: &Vec<Vec2>) -> Result<usize, String> {
    if vertices.len() <= 2 {
        return Err(String::from(
            "vertices vec must contain at least 3 elements",
        ));
    }

    let mut longest_len_squared = (vertices[1] - vertices[0]).length_squared();
    let mut longest_edge: usize = 0;

    for (i, v) in vertices.iter().enumerate().skip(1) {
        let next_vertex = vertices[(i + 1) % vertices.len()]; //
        let len_squared = (next_vertex - *v).length_squared();

        if len_squared < longest_len_squared {
            continue;
        }

        longest_edge = i;
        longest_len_squared = len_squared;
    }

    Ok(longest_edge)
}

pub fn line_intersects_line(a1: Vec2, b1: Vec2, a2: Vec2, b2: Vec2) -> (bool, Vec2) {
    let mut intersection = Vec2::default();

    let x1 = a1.x;
    let x2 = b1.x;
    let x3 = a2.x;
    let x4 = b2.x;
    let y1 = a1.y;
    let y2 = b1.y;
    let y3 = a2.y;
    let y4 = b2.y;

    let d = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
    if d == 0. {
        return (false, intersection);
    }

    let pre = x1 * y2 - y1 * x2;
    let post = x3 * y4 - y3 * x4;
    intersection.x = (pre * (x3 - x4) - (x1 - x2) * post) / d;
    intersection.y = (pre * (y3 - y4) - (y1 - y2) * post) / d;

    (true, intersection)
}

pub fn line_intersects_line_segment(a1: Vec2, b1: Vec2, a2: Vec2, b2: Vec2) -> (bool, Vec2) {
    let (does_intersect, intersection) = line_intersects_line(a1, b1, a2, b2);
    if !does_intersect {
        return (false, intersection);
    }

    let min = Vec2::new(a2.x.min(b2.x), a2.y.min(b2.y));
    let max = Vec2::new(a2.x.max(b2.x), a2.y.max(b2.y));

    (aabb_contains_point(min, max, intersection), intersection)
}

pub fn line_segment_intersects_line_segment(
    a1: Vec2,
    b1: Vec2,
    a2: Vec2,
    b2: Vec2,
) -> (bool, Vec2) {
    let (does_intersect, intersection) = line_intersects_line(a1, b1, a2, b2);
    if !does_intersect {
        return (false, intersection);
    }

    let min = Vec2::new(a1.x.min(b1.x), a1.y.min(b1.y));
    let max = Vec2::new(a1.x.max(b1.x), a1.y.max(b1.y));

    (aabb_contains_point(min, max, intersection), intersection)
}

pub fn get_circle_from_triangle(vertices: &Vec<Vec2>) -> Result<(Vec2, f32), String> {
    if vertices.len() != 3 {
        return Err(String::from("vertices vec must contain exactly 3 elements"));
    }

    let longest_edge = get_longest_edge_of_polygon(vertices).unwrap();
    let a = (longest_edge + 1) % 3;
    let b = (longest_edge + 2) % 3;
    let a_to_b_middle = (vertices[a] + vertices[b]) / 2.;
    let b_to_c_middle = (vertices[b] + vertices[longest_edge]) / 2.;
    let perp_a = (Vec2::new(0., 1.)).rotate(a_to_b_middle - vertices[a]) + a_to_b_middle;
    let perp_b = (Vec2::new(0., 1.)).rotate(b_to_c_middle - vertices[b]) + b_to_c_middle;

    let (_, circle_center) = line_intersects_line(a_to_b_middle, perp_a, b_to_c_middle, perp_b);
    let circle_radius = (vertices[0] - circle_center).length();

    Ok((circle_center, circle_radius))
}

pub fn get_bounding_circle(points: &Vec<Vec2>) -> Result<(Vec2, f32), String> {
    let mut points_on_circle: Vec<Vec2> = Vec::with_capacity(3);
    welzl(points.len(), 0, &points, &mut points_on_circle)
}

fn welzl(
    num_unchecked: usize,
    num_points_on_circle: usize,
    points: &Vec<Vec2>,
    points_on_circle: &mut Vec<Vec2>,
) -> Result<(Vec2, f32), String> {
    if num_unchecked <= 0 || num_points_on_circle == 3 {
        calculate_circle(num_points_on_circle, points_on_circle)
    } else {
        let p = points[num_unchecked - 1];
        match welzl(
            num_unchecked - 1,
            num_points_on_circle,
            &points,
            points_on_circle,
        ) {
            Ok((center, radius)) => {
                if circle_contains_point(center, radius, p) {
                    Ok((center, radius))
                } else {
                    if points_on_circle.len() <= num_points_on_circle {
                        points_on_circle.push(p);
                    } else {
                        points_on_circle[num_points_on_circle] = p;
                    }
                    welzl(
                        num_unchecked - 1,
                        num_points_on_circle + 1,
                        points,
                        points_on_circle,
                    )
                }
            }
            Err(e) => Err(e),
        }
    }
}

fn calculate_circle(
    num_points_on_circle: usize,
    points_on_circle: &Vec<Vec2>,
) -> Result<(Vec2, f32), String> {
    match num_points_on_circle {
        0 => Ok((Vec2::new(0., 0.), 0.)),
        1 => Ok((points_on_circle[0], 0.)),
        2 => {
            let c = (points_on_circle[0] + points_on_circle[1]) / 2.;
            let r = (c - points_on_circle[0]).length();
            Ok((c, r))
        }
        3 => Ok(get_circle_from_triangle(points_on_circle).unwrap()),
        _ => Err(format!(
            "num_points_on_circle is invalid: {}",
            num_points_on_circle
        )),
    }
}

pub fn polygon_contains_point(vertices: &Vec<Vec2>, point: Vec2) -> bool {
    let first_side = point_is_on_left_side_of_line(vertices[0], vertices[1], point);
    for i in 1..vertices.len() {
        let next_index = (i + 1) % vertices.len();
        let side = point_is_on_left_side_of_line(vertices[i], vertices[next_index], point);
        if side != first_side {
            return false;
        }
    }

    true
}

pub fn polygon_intersects_aabb(vertices: &Vec<Vec2>, aabb_min: Vec2, aabb_max: Vec2) -> bool {
    for vertex in vertices {
        if aabb_contains_point(aabb_min, aabb_max, *vertex) {
            return true;
        }
    }

    let aabb_vertices = get_vertices_aabb(aabb_min, aabb_max);
    let aabb_vertices_len = aabb_vertices.len();
    for i in 0..aabb_vertices_len {
        if polygon_contains_point(vertices, aabb_vertices[i]) {
            return true;
        }
    }

    for i in 0..vertices.len() {
        let a = vertices[i];
        let b = vertices[(i + 1) % vertices.len()];
        for j in 0..aabb_vertices_len {
            let other_a = aabb_vertices[j];
            let other_b = aabb_vertices[(j + 1) % vertices.len()];
            if line_segment_intersects_line_segment(a, b, other_a, other_b).0 {
                return true;
            }
        }
    }

    false
}

pub fn polygon_intersects_circle(vertices: &Vec<Vec2>, center: Vec2, radius: f32) -> bool {
    if polygon_contains_point(vertices, center) {
        return true;
    }

    for i in 0..vertices.len() {
        let next_index = (i + 1) % vertices.len();
        let closest_point =
            closest_point_on_line_segment(vertices[i], vertices[next_index], center);
        if circle_contains_point(center, radius, closest_point) {
            return true;
        }
    }

    false
}

pub fn polygon_intersects_polygon(vertices_a: &Vec<Vec2>, vertices_b: &Vec<Vec2>) -> bool {
    for i in 0..vertices_a.len() {
        if polygon_contains_point(vertices_b, vertices_a[i]) {
            return true;
        }
    }

    for i in 0..vertices_b.len() {
        if polygon_contains_point(vertices_a, vertices_b[i]) {
            return true;
        }
    }

    for i in 0..vertices_a.len() {
        let next_index_a = (i + 1) % vertices_a.len();
        for j in 0..vertices_b.len() {
            let next_index_b = (j + 1) % vertices_b.len();
            if line_segment_intersects_line_segment(
                vertices_a[i],
                vertices_a[next_index_a],
                vertices_b[j],
                vertices_b[next_index_b],
            )
            .0
            {
                return true;
            }
        }
    }

    false
}

pub fn polygon_from_aabb(min: Vec2, max: Vec2) -> Vec<Vec2> {
    vec![max, Vec2::new(max.x, min.y), min, Vec2::new(min.x, max.y)]
}

pub fn polygon_from_circle(center: Vec2, radius: f32, num_vertices: usize) -> Vec<Vec2> {
    let mut vertices: Vec<Vec2> = Vec::with_capacity(num_vertices);
    let angle_step = PI * 2. / (num_vertices as f32);
    for i in 0..num_vertices {
        vertices.push(
            center
                + radius * Vec2::new((i as f32 * angle_step).cos(), (i as f32 * angle_step).sin()),
        );
    }

    vertices
}

#[test]
fn calculate_circle_test() {
    // equilateral triangle, this took way longer to figure out than I'd like to admit 🥲
    let points = vec![
        Vec2::new(1., 1.),
        Vec2::new(2., 2.),
        Vec2::new(2.366025, 0.6339746),
    ];

    // calculate_circle should return an error if num_points_on_circle is any value greater than 3
    match calculate_circle(4, &points) {
        Err(_) => assert!(true),
        Ok(_) => assert!(false),
    }

    let (center_a, radius_a) = calculate_circle(1, &points).unwrap();
    assert!(center_a == points[0]);
    assert!(radius_a == 0.);

    let (center_b, radius_b) = calculate_circle(2, &points).unwrap();
    assert!(center_b == Vec2::new(1.5, 1.5));
    assert!(radius_b == (Vec2::new(1.5, 1.5) - points[0]).length());

    let (center_c, radius_c) = calculate_circle(3, &points).unwrap();
    // since this is a equilateral triangle, we can easily calculate the center as the average of all the points
    let expected_center = (points[0] + points[1] + points[2]) / 3.;
    let expected_radius = (expected_center - points[0]).length();
    // there's going to be some rounding errors since these will be calculated differently
    assert!((center_c - expected_center).length() <= 1.68587391E-7);
    assert!(f32::abs(expected_radius - radius_c) <= 1.1920929E-7);
}
