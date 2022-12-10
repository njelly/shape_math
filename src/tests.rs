use super::*;
use bevy::prelude::Vec2;

#[test]
fn aabb_contains_point_test() {
    let min = Vec2::new(0., 0.);
    let max = Vec2::new(1., 1.);

    // all aabb's should contain their own min and max
    assert!(aabb_contains_point(min, max, min));
    assert!(aabb_contains_point(min, max, max));

    // point inside
    assert!(aabb_contains_point(min, max, Vec2::new(0.5, 0.5)));

    // point outside to the left
    assert!(!aabb_contains_point(min, max, Vec2::new(-1., 0.5)));

    // point outside to the right
    assert!(!aabb_contains_point(min, max, Vec2::new(2., 0.5)));

    // point outside above
    assert!(!aabb_contains_point(min, max, Vec2::new(0.5, 2.)));

    // point outside below
    assert!(!aabb_contains_point(min, max, Vec2::new(0.5, -1.)));
}

#[test]
fn aabb_intersects_aabb_test() {
    let min_a = Vec2::new(0., 0.);
    let max_a = Vec2::new(1., 1.);
    let min_b = Vec2::new(1.5, 1.5);
    let max_b = Vec2::new(2., 2.);
    let min_c = Vec2::new(-0.5, -0.5);
    let max_c = Vec2::new(0.5, 0.5);

    // an aabb should intersect itself
    assert!(aabb_intersects_aabb(min_a, max_a, min_a, max_a));

    // these aabbs do not intersect
    assert!(!aabb_intersects_aabb(min_a, max_a, min_b, max_b));
    assert!(!aabb_intersects_aabb(min_c, max_c, min_b, max_b));

    // these aabbs do intersect
    assert!(aabb_intersects_aabb(min_a, max_a, min_c, max_c))
}

#[test]
fn get_vertices_aabb_test() {
    let min = Vec2::new(0., 0.);
    let max = Vec2::new(2.5, 2.5);
    let vertices = get_vertices_aabb(min, max);
    assert!(vertices[0] == max);
    assert!(vertices[1] == Vec2::new(max.x, min.y));
    assert!(vertices[2] == min);
    assert!(vertices[3] == Vec2::new(min.x, max.y));
}

#[test]
fn get_bounding_aabb_test() {
    let min_a = Vec2::new(0., 0.);
    let max_a = Vec2::new(1., 1.);
    let min_b = Vec2::new(1.5, 1.5);
    let max_b = Vec2::new(2., 2.);
    let min_c = Vec2::new(-0.5, -0.5);
    let max_c = Vec2::new(0.5, 0.5);

    let vertices: Vec<Vec2> = get_vertices_aabb(min_a, max_a)
        .into_iter()
        .chain(
            get_vertices_aabb(min_b, max_b)
                .into_iter()
                .chain(get_vertices_aabb(min_c, max_c).into_iter()),
        )
        .collect();

    let (min, max) = get_bounding_aabb(&vertices);
    assert!(min == min_c);
    assert!(max == max_b);

    // the bounding box should contain every point
    for p in vertices {
        assert!(aabb_contains_point(min, max, p));
    }
}

#[test]
fn circle_contains_point_test() {
    let center = Vec2::new(0., 0.);
    let r = 1.;

    // a circle always contains it's center
    assert!(circle_contains_point(center, r, center));

    // a circle contains a point on it's circumference
    assert!(circle_contains_point(center, r, Vec2::new(r, 0.)));

    // a circle does not contain a point barely beyond its radius
    assert!(!circle_contains_point(
        center,
        r,
        Vec2::new(r + f32::EPSILON, 0.)
    ));
}

#[test]
fn circle_intersects_aabb_test() {
    let min_a = Vec2::new(0., 0.);
    let max_a = Vec2::new(1., 1.);
    let min_b = Vec2::new(1., 1.);
    let max_b = Vec2::new(2., 2.);
    let center = Vec2::new(0.5, 1.1);
    let r: f32 = 0.2;

    assert!(circle_intersects_aabb(center, r, min_a, max_a));
    assert!(!circle_intersects_aabb(center, r, min_b, max_b));
}

#[test]
fn circle_intersects_circle_test() {
    let center_a = Vec2::new(0., 0.);
    let r_a: f32 = 0.5;
    let center_b = Vec2::new(r_a, 0.);
    let center_c = Vec2::new(2., 2.);
    let r_c: f32 = 0.4;

    // a circle intersects itself
    assert!(circle_intersects_circle(&center_a, &r_a, &center_a, &r_a));

    // tangent circles count as intersecting
    assert!(circle_intersects_circle(&center_a, &r_a, &center_b, &r_a));

    assert!(!circle_intersects_circle(&center_a, &r_a, &center_c, &r_c))
}

#[test]
fn closest_point_on_line_segment_test() {
    let a = Vec2::new(0., 0.);
    let b = Vec2::new(1., 1.);
    let c = Vec2::new(0.5, 0.5); // a point half-way between a and b
    let d = c + Vec2::new(-1., 1.) * 0.1; // a point a bit off the line perpandicular to ab

    // the points that define a line segment ought to be the closest point on the line segment
    assert!(closest_point_on_line_segment(a, b, a) == a);
    assert!(closest_point_on_line_segment(a, b, b) == b);

    assert!(closest_point_on_line_segment(a, b, c) == c);
    assert!(closest_point_on_line_segment(a, b, d) == c);
    assert!(closest_point_on_line_segment(a, b, Vec2::new(f32::MAX, f32::MAX)) == b);
    assert!(closest_point_on_line_segment(a, b, Vec2::new(f32::MIN, f32::MIN)) == a);
}

#[test]
fn point_is_on_left_side_of_line_test() {
    let a = Vec2::new(0., 0.);
    let b = Vec2::new(0.5, 1.);

    assert!(point_is_on_left_side_of_line(a, b, Vec2::new(0., 1.)));
    assert!(!point_is_on_left_side_of_line(a, b, Vec2::new(1., 0.)));
}

#[test]
fn get_longest_edge_of_polygon_test() {
    // this should fail because the vertices vec has less than 3 elements
    match get_longest_edge_of_polygon(&vec![Vec2::new(0., 0.), Vec2::new(0., 1.)]) {
        Ok(_) => assert!(false),
        Err(_) => assert!(true),
    };

    // the longest edge starts with the point on the second index
    match get_longest_edge_of_polygon(&vec![
        Vec2::new(0., 0.),
        Vec2::new(0.5, 0.2),
        Vec2::new(1., 0.),
    ]) {
        Ok(i) => assert!(i == 2),
        Err(_) => assert!(false),
    }
}

#[test]
fn line_intersects_line_test() {
    let a_start = Vec2::new(0., 0.);
    let a_end = Vec2::new(0., 1.);
    let b_start = Vec2::new(1., 0.);
    let b_end = Vec2::new(1., 3.);
    let c_start = Vec2::new(0.5, 0.5);
    let c_end = Vec2::new(3., 3.);

    let (a_intersects_b, _) = line_intersects_line(a_start, a_end, b_start, b_end);
    let (b_intersects_a, _) = line_intersects_line(b_start, b_end, a_start, a_end);
    let (a_intersects_c, a_intersects_c_point) =
        line_intersects_line(a_start, a_end, c_start, c_end);
    let (c_intersects_a, c_intersects_a_point) =
        line_intersects_line(c_start, c_end, a_start, a_end);
    let (b_intersects_c, b_intersects_c_point) =
        line_intersects_line(b_start, b_end, c_start, c_end);
    let (c_intersects_b, c_intersects_b_point) =
        line_intersects_line(c_start, c_end, b_start, b_end);

    assert!(!a_intersects_b);
    assert!(!b_intersects_a);
    assert!(a_intersects_c);
    assert!(c_intersects_a);
    assert!(a_intersects_c_point == Vec2::new(0., 0.));
    assert!(a_intersects_c_point == c_intersects_a_point);
    assert!(b_intersects_c);
    assert!(c_intersects_b);
    assert!(b_intersects_c_point == Vec2::new(1., 1.));
    assert!(b_intersects_c_point == c_intersects_b_point);
}

#[test]
fn get_circle_from_triangle_test() {
    let triangle = vec![
        Vec2::new(0.0, 0.0),
        Vec2::new(0.5, 0.86602540378),
        Vec2::new(1., 0.),
    ];

    let (circle_center, circle_radius) = get_circle_from_triangle(&triangle).unwrap();
    assert!((circle_center - Vec2::new(0.5, 0.28867513)).length() <= f32::EPSILON);
    assert!((circle_radius - 0.5773502).abs() <= f32::EPSILON);
}

#[test]
fn get_bounding_circle_test() {
    let equilateral_triangle = vec![
        Vec2::new(0.0, 0.0),
        Vec2::new(0.5, 0.86602540378),
        Vec2::new(1., 0.),
    ];

    let some_points = vec![
        Vec2::new(0.0, 0.0),
        Vec2::new(0.5, 0.86602540378),
        Vec2::new(0.5, 0.5), // this one shouldn't affect the bounding circle
        Vec2::new(0.1, 0.4), // ditto
        Vec2::new(1., 0.),
    ];

    let (circle_center_a, circle_radius_a) =
        get_circle_from_triangle(&equilateral_triangle).unwrap();
    let (circle_center_b, circle_radius_b) = get_bounding_circle(&some_points).unwrap();

    // these are calculated differently so there will be rounding errors
    let diff_radius = (circle_radius_a - circle_radius_b).abs();
    let diff_center = (circle_center_a - circle_center_b).length();
    assert!(diff_radius <= 5.96046448E-8);
    assert!(diff_center <= 6.66400197E-8);
}

#[test]
fn polygon_contains_point_test() {
    let triangle = vec![
        Vec2::new(0.0, 0.0),
        Vec2::new(0.5, 0.86602540378),
        Vec2::new(1., 0.),
    ];

    let p_inside = Vec2::new(0.5, 0.5);
    let p_outside = Vec2::new(999., 999.);

    assert!(polygon_contains_point(&triangle, p_inside));
    assert!(!polygon_contains_point(&triangle, p_outside));
}
