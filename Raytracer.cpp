#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iostream>

#include <Vect.h>

#define OBJS_NUM 4

using namespace std;

void writePPM(const char * imagepath, int width, int height, int*** img, int color = 255);

int main()
{
    int width = 256*2, height = 256*2;
    double freq = 3;
    double thickness = 0.2;
    double f = 256;
    double proj_x = 3;
    double e_x = height, e_y = 80, e_z = height;
    Vect eye_pos = Vect(height, 80, height);

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
//    double pn3_x = 0, pn3_y = 0, pn3_z = 1, p3_D = 0;
    Vect pn3 = Vect(0, 0, 1);
    double p3_D = 0;
//    double pn2_x = 1, pn2_y = 0, pn2_z = 0,
    Vect pn2 = Vect(1, 0, 0);
    double p2_D = 0;
//    double pn1_x = 0, pn1_y = 1, pn1_z = 0, p1_D = 0;
    Vect pn1 = Vect(0, 1, 0);
    double p1_D = 0;

//    double sc1_x = 30, sc1_y = 81, sc1_z = 30, s1_r = 50;
    Vect sc1 = Vect(30, 81, 30);
    double s1_r = 50;

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

            double x = (((j - 0)*(f - proj_x - proj_x)) / (height - 1)) + proj_x;
//            double y = height*0.5-i;
            double y = (((i - 0)*(f - proj_x - proj_x)) / (height - 1)) + proj_x;
            double z = -1*(x*n_x + y*n_y + D)/n_z;

            Vect ray_direction = Vect(x - e_x, y - e_y, z - e_z).normalize();

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
                    } else if( object_t_idx == 1 ){
                        intersect_1 = intersect.getVZ();
                        intersect_2 = intersect.getVY()/ height * f;
                    } else if( object_t_idx == 2 ){
                        intersect_1 = intersect.getVX();
                        intersect_2 = intersect.getVY()/ height * f;
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
                        // cout << "aSd" << i << " " << j << endl;
                        p_r = 255;
                        p_g = 0;
                        p_b = 0;
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