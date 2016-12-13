#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <iostream>

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
    int p1_r = 128, p1_g = 0, p1_b = 0;
    int p2_r = 0, p2_g = 128, p2_b = 0;
    int p3_r = 0, p3_g = 0, p3_b = 128;
    int kd = 40, ka = 10;
    // plane normal -- normalized
    double angle = 45;
    double n_x = cos(M_PI * angle / 180), n_y = 0, n_z = sin(M_PI * angle / 180), D = -f * sin(M_PI * angle / 180);
    //double pn3_x = cos(M_PI * angle / 180), pn3_y = 0, pn3_z = -sin(M_PI * angle / 180), p3_D = 35 * sin(M_PI * angle / 180);
    //double pn2_x = -cos(M_PI * angle / 180), pn2_y = 0, pn2_z = -sin(M_PI * angle / 180), p2_D = 35 * sin(M_PI * angle / 180);
    double pn3_x = 0, pn3_y = 0, pn3_z = 1, p3_D = 0;
    double pn2_x = 1, pn2_y = 0, pn2_z = 0, p2_D = 0;
    double pn1_x = 0, pn1_y = 1, pn1_z = 0, p1_D = 0;

    double sc1_x = 30, sc1_y = 81, sc1_z = 30, s1_r = 50;

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

            double dd_x = x - e_x, dd_y = y - e_y;
            double z = -1*(x*n_x + y*n_y + D)/n_z;
            double dd_z = z - e_z;
            double norm = sqrt(dd_x*dd_x + dd_y*dd_y + dd_z*dd_z);
            double d_x = dd_x / norm;
            double d_y = dd_y / norm;
            double d_z = dd_z / norm;

            cout << y << endl;

            double object_ts[OBJS_NUM];
            double object_denoms[OBJS_NUM];

            // planes intersection
            double t1;
            double denom1 = (pn1_x*d_x + pn1_y*d_y + pn1_z*d_z);
            if (abs(denom1) - 1e-5 <= 0)
                t1 = 1e7;
            t1 = -(p1_D + (pn1_x*e_x + pn1_y*e_y + pn1_z*e_z)) / denom1;
            object_ts[0] = t1;
            object_denoms[0] = denom1;
            double t2;
            double denom2 = (pn2_x*d_x + pn2_y*d_y + pn2_z*d_z);
            if (abs(denom2) - 1e-5 <= 0)
                t2 = 1e7;
            t2 = -(p2_D + (pn2_x*e_x + pn2_y*e_y + pn2_z*e_z)) / denom2;
            object_ts[1] = t2;
            object_denoms[1] = denom2;
            double t3;
            double denom3 = (pn3_x*d_x + pn3_y*d_y + pn3_z*d_z);
            if (abs(denom3) - 1e-5 <= 0)
                t3 = 1e7;
            t3 = -(p3_D + (pn3_x*e_x + pn3_y*e_y + pn3_z*e_z)) / denom3;
            object_ts[2] = t3;
            object_denoms[2] = denom3;

            // Sphere intersection
            // Right Sphere
            //
            double C = ( (e_x-sc1_x)*(e_x-sc1_x) + (e_y-sc1_y)*(e_y-sc1_y) + (e_z-sc1_z)*(e_z-sc1_z) ) - s1_r*s1_r;
            double B = 2 * (d_x*(e_x-sc1_x) + d_y*(e_y-sc1_y) + d_z*(e_z-sc1_z));
            double A = d_x*d_x + d_y*d_y + d_z*d_z;
            double D = B*B - 4 * A * C;
            double t4, t4_1, t4_2;
            if ( D < 0 ){
                // No intersection
                t4 = 1e7;
            } else if ( abs(D) < 1e-5 ) {
                // One intersection
                // cout << "1" << endl;
                t4 = (-1*B) / (2*A);
            } else if ( D > 0 ) {
                // Two intersections
                // cout << "2" << endl;
                t4_1 = (-1*B + sqrt(D)) / (2*A);
                t4_2 = (-1*B - sqrt(D)) / (2*A);
                if( t4_1 < t4_2 )
                    t4 = t4_1;
                else
                    t4 = t4_2;
            }
            if( t4 > 0 && t4 < 1e7 ){
                double s_intersect_x = e_x + t4*d_x;
                double s_intersect_y = e_y + t4*d_y;
                double s_intersect_z = e_z + t4*d_z;

                double sn1_x = s_intersect_x - sc1_x, sn1_y = s_intersect_y - sc1_y, sn1_z = s_intersect_z - sc1_z;
                double norm = sqrt(sn1_x*sn1_x + sn1_y*sn1_y + sn1_z*sn1_z);
                sn1_x /= norm;
                sn1_y /= norm;
                sn1_z /= norm;
                double denom4 = (sn1_x*d_x + sn1_y*d_y + sn1_z*d_z);
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
                    double intersect_x = e_x + t*d_x;
                    double intersect_y = e_y + t*d_y;
                    double intersect_z = e_z + t*d_z;
                    double intersect_1, intersect_2;
                    if( object_t_idx == 0 ){
                        intersect_1 = intersect_x;
                        intersect_2 = intersect_z;
                    } else if( object_t_idx == 1 ){
                        intersect_1 = intersect_z;
                        intersect_2 = intersect_y/ height * f;
                    } else if( object_t_idx == 2 ){
                        intersect_1 = intersect_x;
                        intersect_2 = intersect_y/ height * f;
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
