use super::*;
use bevy::prelude::Vec2;

#[test]
fn aabb_contains_point_test() {
    let min = Vec2::new(0., 0.);
    let max = Vec2::new(1., 1.);

    // all aabb's should contain their own min and max
    assert!(aabb_contains_point(&min, &max, &min));
    assert!(aabb_contains_point(&min, &max, &max));

    // point inside
    assert!(aabb_contains_point(&min, &max, &Vec2::new(0.5, 0.5)));

    // point outside to the left
    assert!(!aabb_contains_point(&min, &max, &Vec2::new(-1., 0.5)));

    // point outside to the right
    assert!(!aabb_contains_point(&min, &max, &Vec2::new(2., 0.5)));

    // point outside above
    assert!(!aabb_contains_point(&min, &max, &Vec2::new(0.5, 2.)));

    // point outside below
    assert!(!aabb_contains_point(&min, &max, &Vec2::new(0.5, -1.)));
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
    assert!(aabb_intersects_aabb(&min_a, &max_a, &min_a, &max_a));

    // these aabbs do not intersect
    assert!(!aabb_intersects_aabb(&min_a, &max_a, &min_b, &max_b));
    assert!(!aabb_intersects_aabb(&min_c, &max_c, &min_b, &max_b));

    // these aabbs do intersect
    assert!(aabb_intersects_aabb(&min_a, &max_a, &min_c, &max_c))
}

#[test]
fn get_vertices_aabb_test() {
    let min = Vec2::new(0., 0.);
    let max = Vec2::new(2.5, 2.5);
    let vertices = get_vertices_aabb(&min, &max);
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

    let vertices: Vec<Vec2> = get_vertices_aabb(&min_a, &max_a)
        .into_iter()
        .chain(
            get_vertices_aabb(&min_b, &max_b)
                .into_iter()
                .chain(get_vertices_aabb(&min_c, &max_c).into_iter()),
        )
        .collect();

    let (min, max) = get_bounding_aabb(&vertices);
    assert!(min == min_c);
    assert!(max == max_b);

    // the bounding box should contain every point
    for p in vertices {
        assert!(aabb_contains_point(&min, &max, &p));
    }
}
