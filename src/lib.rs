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

/*public static bool LineIntersectsLine(Vector2 a1, Vector2 b1, Vector2 a2, Vector2 b2,
    out Vector2 intersection, float tolerance = 0.0001f)
{
    intersection = default;

    var x1 = a1.X;
    var x2 = b1.X;
    var x3 = a2.X;
    var x4 = b2.X;

    var y1 = a1.Y;
    var y2 = b1.Y;
    var y3 = a2.Y;
    var y4 = b2.Y;

    // equations of the form x = c (two vertical lines)
    if (MathF.Abs(x1 - x2) < tolerance && MathF.Abs(x3 - x4) < tolerance && MathF.Abs(x1 - x3) < tolerance)
        return false;

    //equations of the form y=c (two horizontal lines)
    if (MathF.Abs(y1 - y2) < tolerance && MathF.Abs(y3 - y4) < tolerance && MathF.Abs(y1 - y3) < tolerance)
        return false;

    //equations of the form x=c (two vertical parallel lines)
    if (MathF.Abs(x1 - x2) < tolerance && MathF.Abs(x3 - x4) < tolerance)
        return false;

    //equations of the form y=c (two horizontal parallel lines)
    if (MathF.Abs(y1 - y2) < tolerance && MathF.Abs(y3 - y4) < tolerance)
        return false;

    //general equation of line is y = mx + c where m is the slope
    //assume equation of line 1 as y1 = m1x1 + c1
    //=> -m1x1 + y1 = c1 ----(1)
    //assume equation of line 2 as y2 = m2x2 + c2
    //=> -m2x2 + y2 = c2 -----(2)
    //if line 1 and 2 intersect then x1=x2=x & y1=y2=y where (x,y) is the intersection p
    //so we will get below two equations
    //-m1x + y = c1 --------(3)
    //-m2x + y = c2 --------(4)

    float x, y;

    //lineA is vertical x1 = x2
    //slope will be infinity
    //so lets derive another solution
    if (MathF.Abs(x1 - x2) < tolerance)
    {
        //compute slope of line 2 (m2) and c2
        var m2 = (y4 - y3) / (x4 - x3);
        var c2 = -m2 * x3 + y3;

        //equation of vertical line is x = c
        //if line 1 and 2 intersect then x1=c1=x
        //subsitute x=x1 in (4) => -m2x1 + y = c2
        // => y = c2 + m2x1
        x = x1;
        y = c2 + m2 * x1;
    }
    //other is vertical x3 = x4
    //slope will be infinity
    //so lets derive another solution
    else if (MathF.Abs(x3 - x4) < tolerance)
    {
        //compute slope of line 1 (m1) and c2
        var m1 = (y2 - y1) / (x2 - x1);
        var c1 = -m1 * x1 + y1;

        //equation of vertical line is x = c
        //if line 1 and 2 intersect then x3=c3=x
        //subsitute x=x3 in (3) => -m1x3 + y = c1
        // => y = c1 + m1x3
        x = x3;
        y = c1 + m1 * x3;
    }
    //lineA & other are not vertical
    //(could be horizontal we can handle it with slope = 0)
    else
    {
        //compute slope of line 1 (m1) and c2
        var m1 = (y2 - y1) / (x2 - x1);
        var c1 = -m1 * x1 + y1;

        //compute slope of line 2 (m2) and c2
        var m2 = (y4 - y3) / (x4 - x3);
        var c2 = -m2 * x3 + y3;

        //solving equations (3) & (4) => x = (c1-c2)/(m2-m1)
        //plugging x value in equation (4) => y = c2 + m2 * x
        x = (c1 - c2) / (m2 - m1);
        y = c2 + m2 * x;

        //verify by plugging intersection p (x, y)
        //in orginal equations (1) & (2) to see if they intersect
        //otherwise x,y values will not be finite and will fail this check
        if (!(MathF.Abs(-m1 * x + y - c1) < tolerance
              && MathF.Abs(-m2 * x + y - c2) < tolerance))
        {
            //return default (no intersection)
            return false;
        }
    }

    //x,y can intersect outside the line segment since line is infinitely long
    //so finally check if x, y is within both the line segments
    intersection = new Vector2(x, y);
    return true;
} */

//pub fn get_circle_from_triangle(points: &Vec<Vec2>) -> Result<(Vec2, f32), String> {
//    if points.len() != 3 {
//        return Err(String::from("vertices vec must contain exactly 3 elements"));
//    }
//
//    let longest_edge = get_longest_edge_of_polygon(points).unwrap();
//    let a = (longest_edge + 1) % 3;
//    let b = (longest_edge + 2) % 3;
//    let a_to_b_middle = (points[a] + points[b]) / 2.;
//    let b_to_c_middle = (points[b] + points[longest_edge]) / 2.;
//    let perp_a = (Vec2::new(0., 1.)).rotate(a_to_b_middle - points[a]) + a_to_b_middle;
//    let perp_b = (Vec2::new(0., 1.)).rotate(b_to_c_middle - points[b]) + b_to_c_middle;
//}

/*
public static unsafe void GetCircleFromTriangleUnsafe(Vector2* points, out Vector2 circleCenter,
    out float circleRadius)
{
    GetLongestEdgeOfPolygonUnsafe(points, 3, out var longestEdge);

    var a = (longestEdge + 1) % 3;
    var b = (longestEdge + 2) % 3;
    var aToBMiddle = (points[a] + points[b]) / 2f;
    var bToCMiddle = (points[b] + points[longestEdge]) / 2f;
    var perpA = (aToBMiddle - points[a]).RotatedByRadians(MathF.PI / 2f) + aToBMiddle;
    var perpB = (bToCMiddle - points[b]).RotatedByRadians(MathF.PI / 2f) + bToCMiddle;

    LineIntersectsLine(aToBMiddle, perpA, bToCMiddle, perpB, out circleCenter);

    circleRadius = (points[0] - circleCenter).Length();
} */

/*
pub fn get_bounding_circle(points: Vec<Vec2>) -> (Vec2, f32) {
    let points_on_cicle: Vec<Vec2> = Vec::with_capacity(3);

    fn welzl(num_unchecked: usize, num_points_on_cicle: usize) -> (Vec2, f32) {
        if num_unchecked == 0 || num_points_on_cicle == 3 {
            return calculate_points_on_circle(num_points_on_cicle);
        }
    }

    fn calculate_circle(num_points_on_circle: usize) -> (Vec2, f32) {
        let c = Vec2::default();
        let r: f32 = 0.;

        return match num_points_on_circle {
            1 => (points_on_circle[0], 0.),
            2 => ((points_on_circle[1] + points_on_circle[0]) / .2, (points_on_circle[1] - points_on_circle[0]).len() / .2),
            3 =>
        }
    }
}
*/

/*
       /// <summary>
       /// Implements Welzl's algorithm for finding the smallest bounding circle containing a set of points in O(n) time.
       /// This website was very helpful: http://www.sunshine2k.de/coding/java/Welzl/Welzl.html
       /// </summary>
       public static unsafe void GetBoundingCircleUnsafe(Vector2* points, int length, out Vector2 circleCenter,
           out float circleRadius)
       {
           var pointsOnCircle = stackalloc Vector2[3];

           welzl(length, 0, out circleCenter, out circleRadius);

           void welzl(int numUnchecked, int numPointsOnCircle, out Vector2 c, out float r)
           {
               if (numUnchecked <= 0 || numPointsOnCircle == 3)
               {
                   calculateCircle(numPointsOnCircle, out c, out r);
                   return;
               }

               var p = points[numUnchecked - 1];
               welzl(numUnchecked - 1, numPointsOnCircle, out c, out r);
               if (CircleContainsPoint(c, r, p))
                   return;

               pointsOnCircle[numPointsOnCircle] = p;
               welzl(numUnchecked - 1, numPointsOnCircle + 1, out c, out r);
           }

           void calculateCircle(int numPointsOnCircle, out Vector2 c, out float r)
           {
               c = default;
               r = default;
               switch (numPointsOnCircle)
               {
                   case 1:
                       c = pointsOnCircle[0];
                       r = 0f;
                       break;
                   case 2:
                       c = (pointsOnCircle[1] + pointsOnCircle[0]) / 2f;
                       r = (pointsOnCircle[1] - pointsOnCircle[0]).Length() / 2f;
                       break;
                   case 3:
                       GetCircleFromTriangleUnsafe(pointsOnCircle, out c, out r);
                       break;
               }
           }
       }
*/
