#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iostream>

#include <Vect.h>

#define OBJS_NUM 6

using namespace std;

void writePPM(const char * imagepath, int width, int height, int*** img, int color = 255);

int main()
{
    int width = 256*2, height = 256*2;
//    double freq = 3;
    double freq = 200;
    double thickness = 0.2;
    double f = 256;
    double proj_x = 3;

    Vect O (0,0,0);
    Vect X (1,0,0);
    Vect Y (0,1,0);
    Vect Z (0,0,1);

//    double e_x = height, e_y = 80, e_z = height;
//    Vect eye_pos = Vect(height, 80, height);
    Vect eye_pos = Vect(5, 0.5, 5);
    Vect look_at = Vect(0, 0.5, 0);
    // Camera three axes are cam_dir, cam_right, cam_down
    Vect cam_dir = (look_at - eye_pos).normalize();
    Vect cam_right = Y.crossProduct(cam_dir).normalize();
    Vect cam_down = cam_dir.crossProduct(cam_right).normalize();
    cam_down = cam_down*(-1);

    int p1_r = 128, p1_g = 0, p1_b = 0;
    int p2_r = 0, p2_g = 128, p2_b = 0;
    int p3_r = 0, p3_g = 0, p3_b = 128;
    int kd = 40, ka = 10;
    // plane normal -- normalized
    double angle = 45;
    double n_x = cos(M_PI * angle / 180), n_y = 0, n_z = sin(M_PI * angle / 180), D = -f * sin(M_PI * angle / 180);
    Vect projection_normal = Vect(n_x, n_y, n_z);
    //double pn3_x = cos(M_PI * angle / 180), pn3_y = 0, pn3_z = -sin(M_PI * angle / 180), p3_D = 35 * sin(M_PI * angle / 180);
    //double pn2_x = -cos(M_PI * angle / 180), pn2_y = 0, pn2_z = -sin(M_PI * angle / 180), p2_D = 35 * sin(M_PI * angle / 180);

    Vect pn1 = Vect(0, 1, 0);
    double p1_D = 0;
    Vect pn2 = Vect(1, 0, 0);
    double p2_D = 0;
    Vect pn3 = Vect(0, 0, 1);
    double p3_D = 0;

    // Right Sphere
    double s1_r = 0.4;
    Vect sc1 = Vect(3.0, s1_r + 0.01, 1.5);
    // Center Sphere
    double s2_r = 0.4;
    Vect sc2 = Vect(2.5, s2_r + 0.01, 2.5);
    // Left Sphere
    double s3_r = 0.4;
    Vect sc3 = Vect(1.5, s3_r + 0.01, 3.0);

    int*** arr = new int**[width];
    for (int i = 0; i < width; i++)
    {
        arr[i] = new int*[height];
        for (int j = 0; j < height; j++)
            arr[i][j] = new int[3]();
    }
    for (int i = 0; i < width; i++)
        for (int j = 0; j < height; j++)
        {

//            double x = (((j - 0)*(f - proj_x - proj_x)) / (height - 1)) + proj_x;
////            double y = height*0.5-i;
//            double y = (((i - 0)*(f - proj_x - proj_x)) / (height - 1)) + proj_x;
//            double z = -1*(x*n_x + y*n_y + D)/n_z;

//            Vect ray_direction = Vect(x - eye_pos.getVX(), y - eye_pos.getVY(), z - eye_pos.getVZ()).normalize();
            double iamnt = (i + 0.5)/width;
            double jamnt = ((height - j) + 0.5)/height;
            Vect ray_direction = (cam_dir + cam_right*(jamnt-0.5) +
                                  cam_down*(iamnt-0.5)).normalize();
//            cout << (cam_down*(jamnt-0.5)).getVX() << " "
//                 << (cam_down*(jamnt-0.5)).getVY() << " "
//                 << (cam_down*(jamnt-0.5)).getVZ() << " " << endl;
            double object_ts[OBJS_NUM];
            double object_denoms[OBJS_NUM];

            // planes intersection
            double t1;

            double denom1 = pn1.dotProduct(ray_direction);
            if (abs(denom1) - 1e-5 <= 0)
                t1 = 1e7;
            t1 = -(p1_D + pn1.dotProduct(eye_pos)) / denom1;
            object_ts[0] = t1;
            object_denoms[0] = denom1;
            double t2;
            double denom2 = pn2.dotProduct(ray_direction);
            if (abs(denom2) - 1e-5 <= 0)
                t2 = 1e7;
            t2 = -(p2_D + pn2.dotProduct(eye_pos)) / denom2;
            object_ts[1] = t2;
            object_denoms[1] = denom2;
            double t3;
            double denom3 = pn3.dotProduct(ray_direction);
            if (abs(denom3) - 1e-5 <= 0)
                t3 = 1e7;
            t3 = -(p3_D + pn3.dotProduct(eye_pos)) / denom3;
            object_ts[2] = t3;
            object_denoms[2] = denom3;

            // Sphere intersection
            // Right Sphere
            double C =  ((eye_pos-sc1).squaredMagnitude()) - s1_r*s1_r;
            double B = ray_direction.dotProduct( eye_pos-sc1 ) * 2;
            double A = ray_direction.squaredMagnitude();
            double D = B*B - 4 * A * C;
            double t4, t4_1, t4_2;
            if ( D < 0 ){
                // No intersection
                t4 = 1e7;
            } else if ( abs(D) < 1e-5 ) {
                // One intersection
                t4 = (-1*B) / (2*A);
            } else if ( D > 0 ) {
                // Two intersections
                t4_1 = (-1*B + sqrt(D)) / (2*A);
                t4_2 = (-1*B - sqrt(D)) / (2*A);
                if( t4_1 < t4_2 )
                    t4 = t4_1;
                else
                    t4 = t4_2;
            }
            if( t4 > 0 && t4 < 1e7 ){
                Vect s_intersect = eye_pos + ray_direction*t4;
                Vect sn1 = (s_intersect - sc1).normalize();
                double denom4 = sn1.dotProduct(ray_direction);
                object_denoms[3] = denom4;
            }
            object_ts[3] = t4;

            // Center Sphere
            C =  ((eye_pos-sc2).squaredMagnitude()) - s2_r*s2_r;
            B = ray_direction.dotProduct( eye_pos-sc2 ) * 2;
            A = ray_direction.squaredMagnitude();
            D = B*B - 4 * A * C;
            double t5, t5_1, t5_2;
            if ( D < 0 ){
                // No intersection
                t5 = 1e7;
            } else if ( abs(D) < 1e-5 ) {
                // One intersection
                t5 = (-1*B) / (2*A);
            } else if ( D > 0 ) {
                // Two intersections
                t5_1 = (-1*B + sqrt(D)) / (2*A);
                t5_2 = (-1*B - sqrt(D)) / (2*A);
                if( t5_1 < t5_2 )
                    t5 = t5_1;
                else
                    t5 = t5_2;
            }
            if( t5 > 0 && t5 < 1e7 ){
                Vect s_intersect = eye_pos + ray_direction*t5;
                Vect sn2 = (s_intersect - sc2).normalize();
                double denom5= sn2.dotProduct(ray_direction);
                object_denoms[4] = denom5;
            }
            object_ts[4] = t5;

            // Left Sphere
            C =  ((eye_pos-sc3).squaredMagnitude()) - s3_r*s3_r;
            B = ray_direction.dotProduct( eye_pos-sc3 ) * 2;
            A = ray_direction.squaredMagnitude();
            D = B*B - 4 * A * C;
            double t6, t6_1, t6_2;
            if ( D < 0 ){
                // No intersection
                t6 = 1e7;
            } else if ( abs(D) < 1e-5 ) {
                // One intersection
                t6 = (-1*B) / (2*A);
            } else if ( D > 0 ) {
                // Two intersections
                t6_1 = (-1*B + sqrt(D)) / (2*A);
                t6_2 = (-1*B - sqrt(D)) / (2*A);
                if( t6_1 < t6_2 )
                    t6 = t6_1;
                else
                    t6 = t6_2;
            }
            if( t6 > 0 && t6 < 1e7 ){
                Vect s_intersect = eye_pos + ray_direction*t6;
                Vect sn3 = (s_intersect - sc3).normalize();
                double denom6= sn3.dotProduct(ray_direction);
                object_denoms[5] = denom6;
            }
            object_ts[5] = t6;

            double denom = 0;
            double t = 0;
            int p_r, p_g, p_b;

            bool intersected_any = false;
            for(int object_t_idx = 0; object_t_idx < OBJS_NUM ; object_t_idx++) {
                bool intersected = true;
                double object_t = object_ts[object_t_idx];
                if( object_t < 0 || object_t >= 1e7){
                    intersected = false;
                    continue;
                }
                for(int other_object_t_idx = 0; other_object_t_idx < OBJS_NUM ; other_object_t_idx++) {
                    double other_object_t = object_ts[other_object_t_idx];
                    if( other_object_t > 0 &&  other_object_t < 1e7 && other_object_t < object_t ){
                        intersected = false;
                        break;
                    }
                }
                if( intersected ){
                    intersected_any = true;
                    denom = object_denoms[object_t_idx];
                    t = object_ts[object_t_idx];
                    Vect intersect = eye_pos + ray_direction*t;
                    double intersect_1, intersect_2;
                    if( object_t_idx == 0 ){
                        intersect_1 = intersect.getVX();
                        intersect_2 = intersect.getVZ();
//                        p_r = 255;
//                        p_g = 0;
//                        p_b = 0;
                    } else if( object_t_idx == 1 ){
                        intersect_1 = intersect.getVZ();
                        intersect_2 = intersect.getVY()/ height * f *2;
//                        p_r = 0;
//                        p_g = 255;
//                        p_b = 0;
                    } else if( object_t_idx == 2 ){
                        intersect_1 = intersect.getVX();
                        intersect_2 = intersect.getVY()/ height * f * 2;
//                        p_r = 0;
//                        p_g = 0;
//                        p_b = 255;
                    }
                    // plane stripes
                    if( object_t_idx < 3 ){
                        if (abs(cos(freq / f*intersect_1 * 2 * M_PI)) < thickness || abs(cos(freq / f*intersect_2 * 2 * M_PI)) < thickness)
                        {
                            p_r = 0;
                            p_g = 0;
                            p_b = 0;
                        }
                        else
                        {
                            p_r = 255;
                            p_g = 255;
                            p_b = 255;
                        }
                    } else if( object_t_idx == 3 ) {
                        // Right Sphere
                        p_r = 255;
                        p_g = 0;
                        p_b = 0;
                    } else if( object_t_idx == 4 ) {
                        // Center Sphere
                        p_r = 0;
                        p_g = 255;
                        p_b = 0;
                    } else if( object_t_idx == 5 ) {
                        // Left Sphere
                        p_r = 0;
                        p_g = 0;
                        p_b = 255;
                    }

                    break;
                }
            }

            if( !intersected_any ){
                continue;
            }
            arr[i][j][0] = ka + kd*abs(denom)*p_r;
            arr[i][j][1] = ka + kd*abs(denom)*p_g;
            arr[i][j][2] = ka + kd*abs(denom)*p_b;

            // clip values
            arr[i][j][0] = arr[i][j][0] > 255 ? 255 : arr[i][j][0];
            arr[i][j][1] = arr[i][j][1] > 255 ? 255 : arr[i][j][1];
            arr[i][j][2] = arr[i][j][2] > 255 ? 255 : arr[i][j][2];
        }
    writePPM("test.ppm", width, height, arr);
    return 0;
}

void writePPM(const char * imagepath, int width, int height, int*** img, int color)
{
    std::ofstream f(imagepath);
    // header
    f << "P3" << "\n";
    // meta
    f << width << " " << height << " " << color << "\n";
    // data
    int r = 255, g = 0, b = 0;
    for (int i = 0; i < width; i++)
        for (int j = 0; j < height; j++)
            f << img[i][j][0] << " " << img[i][j][1] << " " << img[i][j][2] << "\n";

    f.close();
}
