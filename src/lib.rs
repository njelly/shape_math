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
