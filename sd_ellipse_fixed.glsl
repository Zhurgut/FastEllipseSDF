




// // create the UBO to upload to GPU

// struct EllipseUBO {
//     float radius[2];
//     float r;
//     float inv_r;
//     float inv_ry;
//     float parabola_peak_x;
//     float minus_inv2A;
//     float A, B, C;
// };


// struct EllipseUBO ellipse_ubo(float rx, float ry) {

//     // assert rx >= ry

//     float r = rx / ry;
//     float inv_r = 1.0f / r;
//     float peak_x = r - inv_r;
//     float A = -r / (peak_x + 1e-6);

//     struct EllipseUBO data;
//     data.radius[0] = rx;
//     data.radius[1] = ry;
//     data.r = r;
//     data.inv_r = inv_r;
//     data.parabola_peak_x = peak_x;
//     data.minus_inv2A = -0.5 * peak_x * inv_r;
//     data.A = A;
//     data.B = 2.0 * r;
//     data.C = 1.0 - (r * r);

//     return data;
// }





layout(std140, binding = 0) uniform EllipseBlock {
    vec2  radius;
    float r;
    float inv_r;
    float inv_ry;
    float parabola_peak_x;
    float minus_inv2A;
    float A, B, C;
} e;


// float fma(float a, float b, float c) {
//     return a*b + c;
// }

float parabola(float A, float B, float C, float x) {
    return fma(fma(A, x, B), x, C);
}


float sd_ellipse(vec2 p) {            
    p = abs(p);
    p *= e.inv_ry;

    float A = e.A;
    float B = e.B;
    float C = e.C;
    
    // float a = -A;
    float b = 2.0*A*p.x;
    float c = fma(B, p.x, C - p.y);

    float disc = max(0.0, b*b + 4.0*A*c);

    float tx = (b + sqrt(disc)) * e.minus_inv2A;
    tx = clamp(tx, 0.0, e.parabola_peak_x - 1e-4);
    float ty = parabola(A, B, C, tx);

    vec2 t = vec2(tx, ty);

    vec2 p_shift = vec2(p.x * e.inv_r, p.y);
    vec2 t_shift = vec2(t.x * e.inv_r, t.y);

    vec2 dif = t_shift - p_shift;

    float minus_2a2 = -2.0*dot(dif, dif);
    float b2 = 2.0 * dot(dif, p_shift);
    float c2 = dot(p_shift, p_shift) - 1.0;

    float disc2 = max(0.0, b2*b2 + 2.0*minus_2a2*c2);

    float z = (b2 + sqrt(disc2)) / (minus_2a2 - 1e-20); 

    vec2 intersect = mix(p_shift, t_shift, z);
    intersect.x = intersect.x * e.r;

    return sign(c2) * e.radius.y * distance(intersect, p);

}



