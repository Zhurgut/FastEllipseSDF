

// float fma(float a, float b, float c) {
//     return a*b + c;
// }

float parabola(float A, float B, float C, float x) {
    return fma(fma(A, x, B), x, C);
}


float sd_ellipse( vec2 p, vec2 radius )
{
    if (radius.x < radius.y ) { 
        p=p.yx; radius=radius.yx; 
    }
    
    float inv_ry = 1.0 / radius.y;

    float r = radius.x * inv_ry;
    float inv_r = 1.0 / r;
    float parabola_peak_x = r - inv_r;
    float minus_inv2A = -0.5 * parabola_peak_x * inv_r;

    float A = -r / (parabola_peak_x + 1e-6);
    float B = 2.0*r;
    float C = fma(-r, r, 1.0);

            
    p = abs(p);
    p *= inv_ry;
    
    // float a = -A;
    float b = 2.0*A*p.x;
    float c = fma(B, p.x, C - p.y);

    float disc = max(0.0, b*b + 4.0*A*c);

    float tx = (b + sqrt(disc)) * minus_inv2A;
    tx = clamp(tx, 0.0, parabola_peak_x - 1e-4);
    float ty = parabola(A, B, C, tx);

    vec2 t = vec2(tx, ty);

    vec2 p_shift = vec2(p.x * inv_r, p.y);
    vec2 t_shift = vec2(t.x * inv_r, t.y);

    vec2 dif = t_shift - p_shift;

    float minus_2a2 = -2.0*dot(dif, dif);
    float b2 = 2.0 * dot(dif, p_shift);
    float c2 = dot(p_shift, p_shift) - 1.0;

    float disc2 = max(0.0, b2*b2 + 2.0*minus_2a2*c2);

    float z = (b2 + sqrt(disc2)) / (minus_2a2 - 1e-20); 

    vec2 intersect = mix(p_shift, t_shift, z);
    intersect.x = intersect.x * r;

    return sign(c2) * radius.y * distance(intersect, p);

}
