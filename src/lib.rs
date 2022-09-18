use bevy::math::Vec2;

#[cfg(test)]
mod tests;

pub fn aabb_contains_point(min: &Vec2, max: &Vec2, p: &Vec2) -> bool {
    p.x >= min.x && p.x <= max.x && p.y >= min.y && p.y <= max.y
}

pub fn aabb_intersects_aabb(min_a: &Vec2, max_a: &Vec2, min_b: &Vec2, max_b: &Vec2) -> bool {
    max_a.x - min_a.x + max_b.x - min_b.x >= max_a.x.max(max_b.x) - min_a.x.min(min_b.x)
        && max_a.y - min_a.y + max_b.y - min_b.y >= max_a.y.max(max_b.y) - min_a.y.min(min_b.y)
}

pub fn get_vertices_aabb(min: &Vec2, max: &Vec2) -> Vec<Vec2> {
    vec![
        max.clone(),
        Vec2::new(max.x, min.y),
        min.clone(),
        Vec2::new(min.x, max.y),
    ]
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

pub fn circle_contains_point(center: &Vec2, radius: &f32, p: &Vec2) -> bool {
    (*p - *center).length_squared() <= radius * radius
}

pub fn circle_intersects_aabb(center: &Vec2, r: &f32, aabb_min: &Vec2, aabb_max: &Vec2) -> bool {
    if aabb_contains_point(&aabb_min, &aabb_max, &center) {
        return true;
    }

    let vertices = get_vertices_aabb(&aabb_min, &aabb_max);
    for (i, v) in vertices.iter().enumerate() {
        let next_index = (i + 1) % vertices.len();
        let closest_point_on_edge =
            closest_point_on_line_segment(v, &vertices[next_index], &center);
        if circle_contains_point(&center, &r, &closest_point_on_edge) {
            return true;
        }
    }

    false
}

pub fn circle_intersects_circle(center_a: &Vec2, r_a: &f32, center_b: &Vec2, r_b: &f32) -> bool {
    let r_sum = r_a + r_b;
    (*center_b - *center_a).length_squared() <= r_sum * r_sum
}

pub fn closest_point_on_line_segment(a: &Vec2, b: &Vec2, p: &Vec2) -> Vec2 {
    let ap = *p - *a;
    let ab = *b - *a;
    let dist = ap.dot(ab) / ab.length_squared();
    if dist <= 0. {
        *a
    } else if dist >= 1. {
        *b
    } else {
        return *a + ab * dist;
    }
}

pub fn point_is_on_left_side_of_line(a: &Vec2, b: &Vec2, p: &Vec2) -> bool {
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

pub fn line_intersects_line(a1: &Vec2, b1: &Vec2, a2: &Vec2, b2: &Vec2) -> (bool, Vec2) {
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

pub fn line_intersects_line_segment(a1: &Vec2, b1: &Vec2, a2: &Vec2, b2: &Vec2) -> (bool, Vec2) {
    let (does_intersect, intersection) = line_intersects_line(a1, b1, a2, b2);
    if !does_intersect {
        return (false, intersection);
    }

    let min = Vec2::new(a2.x.min(b2.x), a2.y.min(b2.y));
    let max = Vec2::new(a2.x.max(b2.x), a2.y.max(b2.y));

    (aabb_contains_point(&min, &max, &intersection), intersection)
}

pub fn line_segment_intersects_line_segment(
    a1: &Vec2,
    b1: &Vec2,
    a2: &Vec2,
    b2: &Vec2,
) -> (bool, Vec2) {
    let (does_intersect, intersection) = line_intersects_line(a1, b1, a2, b2);
    if !does_intersect {
        return (false, intersection);
    }

    let min = Vec2::new(a1.x.min(b1.x), a1.y.min(b1.y));
    let max = Vec2::new(a1.x.max(b1.x), a1.y.max(b1.y));

    (aabb_contains_point(&min, &max, &intersection), intersection)
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

    let (_, circle_center) = line_intersects_line(&a_to_b_middle, &perp_a, &b_to_c_middle, &perp_b);
    let circle_radius = (vertices[0] - circle_center).length();

    Ok((circle_center, circle_radius))
}

pub fn get_bounding_circle(points: &Vec<Vec2>) -> (Vec2, f32) {
    let mut points_on_circle: Vec<Vec2> = Vec::with_capacity(3);

    fn welzl(
        num_unchecked: usize,
        num_points_on_cicle: usize,
        points_on_circle: &mut Vec<Vec2>,
    ) -> (Vec2, f32) {
        if num_unchecked == 0 || num_points_on_cicle == 3 {
            return calculate_circle(num_points_on_cicle, points_on_circle);
        }

        let p = points_on_circle[num_unchecked - 1];
        let (c, r) = welzl(num_unchecked - 1, num_points_on_cicle, points_on_circle);
        if circle_contains_point(&c, &r, &p) {
            return (c, r);
        }

        points_on_circle[num_points_on_cicle] = p;
        welzl(num_unchecked, num_points_on_cicle, points_on_circle)
    }

    fn calculate_circle(num_points_on_circle: usize, points_on_circle: &Vec<Vec2>) -> (Vec2, f32) {
        return match num_points_on_circle {
            1 => (points_on_circle[0], 0.),
            2 => (
                (points_on_circle[1] + points_on_circle[0]) / 2.,
                (points_on_circle[1] - points_on_circle[0]).length() / 2.,
            ),
            _ => get_circle_from_triangle(points_on_circle).unwrap(),
        };
    }

    welzl(points.len(), 0, &mut points_on_circle)
}
