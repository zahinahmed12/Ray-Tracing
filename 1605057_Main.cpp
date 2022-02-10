#include<bits/stdc++.h>
#include "bitmap_image.hpp"
#include "1605057_Header.h"
#include <windows.h>
#include <glut.h>
#define pi (2*acos(0.0))
using namespace  std;

int draw_axes;
int recursion_level;
int pixels;
int no_of_objects;
int no_of_lights;
vector<Object*> objects;
vector<Light*> lights;

struct vector3D cam_pos, u, r, l ;
double move_cam;
double rotate_cam;
const int window_height=500, window_width=500, view_angle=80;

struct vector3D new_vector(struct vector3D a, struct vector3D temp, double angle)
{
    struct vector3D res{};
    double alpha = cos(angle * (pi / 180));
    double beta = sin(angle * (pi / 180));
    res.set_point(alpha * a.x + beta * temp.x,alpha * a.y + beta * temp.y,alpha * a.z + beta * temp.z);
    return res;
}
vector<double> split_line(const string& str)
{
    string word;
    istringstream ss(str);
    vector<double> aux;
    while(ss >> word)
    {
        aux.push_back(stod(word));
    }
    return aux;
}
void drawAxes()
{
    if(draw_axes==1)
    {
        //glColor3f(1.0, 1.0, 1.0);
        glBegin(GL_LINES);{
            glColor3f(1.0, 0, 0);
            glVertex3f( 1000,0,0);
            glVertex3f(-1000,0,0);

            glColor3f(0, 1.0, 0);
            glVertex3f(0,-1000,0);
            glVertex3f(0, 1000,0);

            glColor3f(0, 0, 1.0);
            glVertex3f(0,0, 1000);
            glVertex3f(0,0,-1000);
        }glEnd();
    }
}
void drawSS()
{
    drawAxes();
    for(int i=0; i<no_of_objects; i++)
    {
        glPushMatrix();
        glTranslatef(objects[i]->ref_point.x,objects[i]->ref_point.y,objects[i]->ref_point.z);
        objects[i]->draw();
        glPopMatrix();
    }

    objects[no_of_objects]->draw();

    for(int i=0; i<no_of_lights; i++)
    {
        glTranslatef(lights[i]->light_pos.x,lights[i]->light_pos.y,lights[i]->light_pos.z);
        lights[i]->draw();
    }

}
void loadData()
{
    vector<string> lines;
    string line;

    ifstream scene_file("..\\scene.txt");
    while(getline(scene_file, line)) {
        lines.push_back(line);
    }
    scene_file.close();

    unsigned int no_of_lines = lines.size();
    //cout<<no_of_lines<<endl;

    vector<double> temp1 =  split_line(lines[0]);
    recursion_level = (int)temp1[0];
    temp1.clear();

    temp1 = split_line(lines[1]);
    pixels = (int) temp1[0];
    temp1.clear();

    temp1= split_line(lines[3]);
    no_of_objects = (int) temp1[0];
    temp1.clear();

    int i;
    for(i=4; i<no_of_lines; i++)
    {
        if(lines[i]=="sphere")
        {
            vector<double> aux = split_line(lines[i+1]);
            vector3D p{};
            p.set_point(aux[0],aux[1],aux[2]);
            //cout<<aux[0]<<" "<<aux[1]<<" "<<aux[2]<<" "<<endl;
            vector<double> aux2 = split_line(lines[i+2]);
            //cout<<aux2[0]<<endl;
            Object *sphere = new Sphere(p,aux2[0]);

            aux.clear();
            aux = split_line(lines[i+3]);
            sphere->setColor(aux[0],aux[1],aux[2]);
            //cout<<aux[0]<<" "<<aux[1]<<" "<<aux[2]<<" "<<endl;

            aux.clear();
            aux = split_line(lines[i+4]);
            sphere->setCoefficient(aux[0],aux[1],aux[2],aux[3]);
            //cout<<aux[0]<<" "<<aux[1]<<" "<<aux[2]<<" "<<aux[3]<<endl;

            aux.clear();
            aux = split_line(lines[i+5]);
            sphere->setShine((int)aux[0]);
            //cout<<aux[0]<<endl;

            sphere->obj_type="sphere";
            objects.push_back(sphere);
            i+=6;

            aux.clear();
            aux2.clear();
        }
        else if(lines[i]=="triangle")
        {
            vector<vector3D> points;
            vector3D p{};

            vector<double> aux = split_line(lines[i+1]);
            p.set_point(aux[0],aux[1],aux[2]);
            points.push_back(p);
            aux.clear();

            aux = split_line(lines[i+2]);
            p.set_point(aux[0],aux[1],aux[2]);
            points.push_back(p);
            aux.clear();

            aux = split_line(lines[i+3]);
            p.set_point(aux[0],aux[1],aux[2]);
            points.push_back(p);
            aux.clear();

            Object *triangle = new Triangle(points);

            aux = split_line(lines[i+4]);
            triangle->setColor(aux[0],aux[1],aux[2]);
            aux.clear();

            aux = split_line(lines[i+5]);
            triangle->setCoefficient(aux[0],aux[1],aux[2],aux[3]);
            aux.clear();

            aux = split_line(lines[i+6]);
            triangle->setShine((int)aux[0]);

            triangle->obj_type="triangle";
            objects.push_back(triangle);

            i+=7;
            aux.clear();
        }
        else if(lines[i]=="general")
        {
            vector<double> aux = split_line(lines[i+1]);

            vector3D p{};
            vector<double> aux2 = split_line(lines[i+2]);
            p.set_point(aux2[0],aux2[1],aux2[2]);

            Object *quad = new Quadric(p, aux2[3], aux2[4], aux2[5]);
            quad->set_eqn_val(aux);

            aux2.clear();
            aux2 = split_line(lines[i+3]);
            quad->setColor(aux2[0],aux2[1],aux2[2]);

            aux2.clear();
            aux2 = split_line(lines[i+4]);
            quad->setCoefficient(aux2[0],aux2[1],aux2[2],aux2[3]);

            aux2.clear();
            aux2 = split_line(lines[i+5]);
            quad->setShine((int)aux2[0]);

            quad->obj_type="general";
            objects.push_back(quad);
            i+=6;
            aux.clear();
            aux2.clear();
        }
        else {
            break;
        }
    }
    //cout<<i<<endl;
    vector<double> aux= split_line(lines[i]);
    no_of_lights = (int)aux[0];
    aux.clear();

    for(int j=1;j<no_of_lights*2;j=j+2)
    {
        vector3D p{};
        vector<double> aux2= split_line(lines[i+j]);
        //cout<<aux2[0]<<" "<<aux2[1]<<" "<<aux2[2]<<" "<<endl;
        p.set_point(aux2[0],aux2[1],aux2[2]);
        vector<double> aux3= split_line(lines[i+j+1]);
        //cout<<aux3[0]<<" "<<aux3[1]<<" "<<aux3[2]<<" "<<endl;
        auto *light= new Light(p, aux3[0], aux3[1], aux3[2]);
        lights.push_back(light);
        aux2.clear();
        aux3.clear();
    }

    Object *floor = new Floor(1000, 20);
    floor->obj_type="floor";
    floor->setColor(1,1,1);
    floor->setCoefficient(0.4,0.1,0.3,0.7);
    floor->setShine(10);
    objects.push_back(floor);

    lines.clear();
}
void unLoadData()
{
    for(auto obj:objects)
    {
        delete obj;
    }
    for(auto lit:lights)
    {
        delete lit;
    }
    objects.clear();
    lights.clear();
}
void printObjects()
{
    for(int i=0; i<no_of_objects+1; i++)
    {
        if(objects[i]->obj_type=="sphere")
        {
            cout<<"Sphere\n";
            cout<< "ref_point: "<<objects[i]->ref_point.x <<" "<< objects[i]->ref_point.y <<" "<< objects[i]->ref_point.z<<endl;
            cout<< "radius: "<<objects[i]->length<<endl;
            cout<< "color: "<<objects[i]->color[0] <<" "<< objects[i]->color[1] <<" "<< objects[i]->color[2]<<endl;
            cout<< "coefficient: "<<objects[i]->coEfficients[0] <<" "<< objects[i]->coEfficients[1] <<" "<< objects[i]->coEfficients[2]<<" "<< objects[i]->coEfficients[3]<<endl;
            cout<< "shine: "<<objects[i]->shine<<endl;
        }
        else if(objects[i]->obj_type=="triangle")
        {
            cout<<"Triangle\n";
            cout<< "1st_point: "<<objects[i]->tri_points[0].x <<" "<< objects[i]->tri_points[0].y <<" "<< objects[i]->tri_points[0].z<<endl;
            cout<< "2nd_point: "<<objects[i]->tri_points[1].x <<" "<< objects[i]->tri_points[1].y <<" "<< objects[i]->tri_points[1].z<<endl;
            cout<< "3rd_point: "<<objects[i]->tri_points[2].x <<" "<< objects[i]->tri_points[2].y <<" "<< objects[i]->tri_points[2].z<<endl;
            cout<< "color: "<<objects[i]->color[0] <<" "<< objects[i]->color[1] <<" "<< objects[i]->color[2]<<endl;
            cout<< "coefficient: "<<objects[i]->coEfficients[0] <<" "<< objects[i]->coEfficients[1] <<" "<< objects[i]->coEfficients[2]<<" "<< objects[i]->coEfficients[3]<<endl;
            cout<< "shine: "<<objects[i]->shine<<endl;
        }
        else if(objects[i]->obj_type=="general")
        {
            cout<<"General\n";
            cout<< "ref_point: "<<objects[i]->ref_point.x <<" "<< objects[i]->ref_point.y <<" "<< objects[i]->ref_point.z<<endl;
            cout<< "length: "<<objects[i]->length<<endl;
            cout<< "width: "<<objects[i]->width<<endl;
            cout<< "height: "<<objects[i]->height<<endl;
            cout<< "ABCDEFGHIJ: "<<objects[i]->eqn_val[0] <<" "<< objects[i]->eqn_val[1]<<" "<< objects[i]->eqn_val[2]<<" "<< objects[i]->eqn_val[3]<<" "<< objects[i]->eqn_val[4]<<" "<< objects[i]->eqn_val[5]<<" "<< objects[i]->eqn_val[6]<<" "<< objects[i]->eqn_val[7]<<" "<< objects[i]->eqn_val[8]<<" "<< objects[i]->eqn_val[9]<<endl;
            cout<< "color: "<<objects[i]->color[0] <<" "<< objects[i]->color[1] <<" "<< objects[i]->color[2]<<endl;
            cout<< "coefficient: "<<objects[i]->coEfficients[0] <<" "<< objects[i]->coEfficients[1] <<" "<< objects[i]->coEfficients[2]<<" "<< objects[i]->coEfficients[3]<<endl;
            cout<< "shine: "<<objects[i]->shine<<endl;
        }
        else if(objects[i]->obj_type=="floor")
        {
            cout<<"Floor\n";
            cout<< "ref_point: "<<objects[i]->ref_point.x <<" "<< objects[i]->ref_point.y <<" "<< objects[i]->ref_point.z<<endl;
            cout<< "tile_width: "<<objects[i]->length<<endl;
            //cout<< "color: "<<objects[i]->color[0] <<" "<< objects[i]->color[1] <<" "<< objects[i]->color[2]<<endl;
            //cout<< "coefficient: "<<objects[i]->coEfficients[0] <<" "<< objects[i]->coEfficients[1] <<" "<< objects[i]->coEfficients[2]<<" "<< objects[i]->coEfficients[3]<<endl;
            //cout<< "shine: "<<objects[i]->shine<<endl;
        }
    }
}
void printLights()
{
    for(int i=0;i<no_of_lights;i++)
    {
        cout<< "Light position: "<<lights[i]->light_pos.x <<" "<< lights[i]->light_pos.y <<" "<< lights[i]->light_pos.z<<endl;
        cout<< "color: "<<lights[i]->color[0] <<" "<< lights[i]->color[1] <<" "<< lights[i]->color[2]<<endl;
    }
}

void capture()
{
    bitmap_image image(pixels, pixels);
    image.clear();
    double plane_distance = (window_height/2.0)/tan(view_angle/2.0*(pi/180.0));
    double du = 1.0*window_width/pixels;
    double dv = 1.0*window_height/pixels;
    vector3D top_left{};
    top_left.set_point(cam_pos.x + l.x*plane_distance - r.x*window_width/2 + u.x*window_height/2,
                       cam_pos.y + l.y*plane_distance - r.y*window_width/2 + u.y*window_height/2,
                       cam_pos.z + l.z*plane_distance - r.z*window_width/2 + u.z*window_height/2);

    top_left.set_point(top_left.x + r.x*0.5*du - u.x*0.5*dv,
                       top_left.y + r.y*0.5*du - u.y*0.5*dv,
                       top_left.z + r.z*0.5*du - u.z*0.5*dv);
    //top_left.print_point();

    for(int i=1; i<=pixels; i++)
    {
        for(int j=1; j<=pixels; j++)
        {
            vector3D cur_pixel{}, ray_dir{};
            cur_pixel.set_point(top_left.x + r.x*i*du - u.x*j*dv,
                                top_left.y + r.y*i*du - u.y*j*dv,
                                top_left.z + r.z*i*du - u.z*j*dv);

            ray_dir.set_point(cur_pixel.x-cam_pos.x,cur_pixel.y-cam_pos.y,cur_pixel.z-cam_pos.z);
            Ray *ray = new Ray(cam_pos,ray_dir);

            auto *color = new double [3];
            int nearest=-1;
            double t=max_value, t_min;

            for(int k=0; k<no_of_objects+1; k++)
            {
                double temp = objects[k]->intersect(ray,color,0);

                if(0<temp && temp<t){
                    t = temp;
                    nearest = k ;
                }
            }
            if(nearest!=-1){
                t_min = objects[nearest]->intersect(ray,color,1);
                image.set_pixel(i-1,j-1,(int)(color[0]*255),(int)(color[1]*255),(int)(color[2]*255));
            }
            delete ray;
            delete[] color;
        }
    }
    image.save_image("output.bmp");
}

void keyboardListener(unsigned char key, int x, int y){
    struct vector3D temp1{};
    struct vector3D temp2{};
    switch(key){
        case '1':
            temp1 = cross_product(u,l);
            l = new_vector(l,temp1, rotate_cam);
            temp2 = cross_product(u,r);
            r = new_vector(r,temp2, rotate_cam);

            break;
        case '2':
            temp1 = cross_product(l,u);
            l = new_vector(l,temp1, rotate_cam);
            temp2 = cross_product(r,u);
            r = new_vector(r,temp2,rotate_cam);

            break;
        case '3':
            temp1 = cross_product(r,l);
            l = new_vector(l,temp1,rotate_cam);
            temp2 = cross_product(r,u);
            u = new_vector(u,temp2,rotate_cam);

            break;
        case '4':
            temp1 = cross_product(l,r);
            l = new_vector(l,temp1,rotate_cam);
            temp2 = cross_product(u,r);
            u = new_vector(u,temp2,rotate_cam);

            break;
        case '5':
            temp1 = cross_product(u,l);
            u = new_vector(u,temp1,rotate_cam);
            temp2 = cross_product(r,l);
            r = new_vector(r,temp2,rotate_cam);

            break;
        case '6':
            temp1 = cross_product(l,u);
            u = new_vector(u,temp1,rotate_cam);
            temp2 = cross_product(l,r);
            r = new_vector(r,temp2,rotate_cam);
            break;

        case '0':
            cout<<"Capturing...\n";
            capture();
            cout<<"Done.\n";
            break;

        default:
            break;
    }
    glutPostRedisplay();
}
void specialKeyListener(int key, int x, int y){
    switch(key){
        case GLUT_KEY_DOWN:		//down arrow key
            cam_pos.set_point(cam_pos.x-move_cam * l.x, cam_pos.y-move_cam * l.y, cam_pos.z-move_cam * l.z);
            break;
        case GLUT_KEY_UP:		// up arrow key
            cam_pos.set_point(cam_pos.x+move_cam * l.x, cam_pos.y+move_cam * l.y, cam_pos.z+move_cam * l.z);
            break;
        case GLUT_KEY_RIGHT:
            cam_pos.set_point(cam_pos.x+move_cam * r.x, cam_pos.y+move_cam * r.y, cam_pos.z+move_cam * r.z);
            break;
        case GLUT_KEY_LEFT:
            cam_pos.set_point(cam_pos.x-move_cam * r.x, cam_pos.y-move_cam * r.y, cam_pos.z-move_cam * r.z);
            break;
        case GLUT_KEY_PAGE_UP:
            cam_pos.set_point(cam_pos.x+move_cam * u.x, cam_pos.y+move_cam * u.y, cam_pos.z+move_cam * u.z);
            break;
        case GLUT_KEY_PAGE_DOWN:
            cam_pos.set_point(cam_pos.x-move_cam * u.x, cam_pos.y-move_cam * u.y, cam_pos.z-move_cam * u.z);
            break;
        default:
            break;
    }
    glutPostRedisplay();
}
void display(){

    //clear the display
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0,0,0,0);	//color black
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    /********************
    / set-up camera here
    ********************/
    //load the correct matrix -- MODEL-VIEW matrix
    glMatrixMode(GL_MODELVIEW);

    //initialize the matrix
    glLoadIdentity();

    //now give three info
    //1. where is the camera (viewer)?
    //2. where is the camera looking?
    //3. Which direction is the camera's UP direction?

    //gluLookAt(100,100,100,	0,0,0,	0,0,1);
    //gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
    //gluLookAt(0,0,200,	0,0,0,	0,1,0);
    gluLookAt(cam_pos.x,cam_pos.y,cam_pos.z,	cam_pos.x+l.x,cam_pos.y+l.y,cam_pos.z+l.z,	u.x,u.y,u.z);


    //again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);


    /****************************
    / Add your objects from here
    ****************************/

    drawSS();

    //ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
    glutSwapBuffers();
}
void animate(){
    //glutPostRedisplay();
}
void init(){
    //codes for initialization
    draw_axes=1;

    cam_pos.x = 100;
    cam_pos.y = 100;
    cam_pos.z = 20;
    u.x = 0;
    u.y = 0;
    u.z = 1;
    r.x = -1.0/sqrt(2);
    r.y = 1.0/sqrt(2);
    r.z = 0;
    l.x = -1.0/sqrt(2);
    l.y = -1.0/sqrt(2);
    l.z = 0;

    move_cam=5.0;
    rotate_cam=3.0;

    //clear the screen
    glClearColor(0,0,0,0);

    /************************
    / set-up projection here
    ************************/
    //load the PROJECTION matrix
    glMatrixMode(GL_PROJECTION);

    //initialize the matrix
    glLoadIdentity();

    //give PERSPECTIVE parameters
    gluPerspective(view_angle,	1,	1,	1000.0);
    //field of view in the Y (vertically)
    //aspect ratio that determines the field of view in the X direction (horizontally)
    //near distance
    //far distance
}

int main(int argc, char **argv){

    loadData();
    atexit(unLoadData);
//    printObjects();
//    printLights();
//    cout<<recursion_level<<" "<<pixels<<" "<<no_of_objects<<" "<<no_of_lights<<endl;
//    cout<<objects.size()<<endl;
//    check();

    glutInit(&argc,argv);
    glutInitWindowSize(window_width, window_height);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

    glutCreateWindow("Ray Tracing OpenGL Program 1");

    init();
    //capture();

//    vector3D center{};
//    center.set_point(10,4,-5);
//    Object *sphere = new Sphere(center,10.0);
//    sphere->setColor(1,0,0);
//    sphere->setCoefficient(0.4,0.2,0.2,0.2);
//    vector3D start{}, dir{};
//    start.set_point(10,0,0);
//    dir.set_point(10,3,-1);
//    Ray *rr = new Ray(start,dir);
//    rr->normalized_dir();
//    auto *color = new double [3];
//    sphere->intersect(rr,color,1);
//    cout<<color[0]<<" "<<color[1]<<" "<<color[2];

    glEnable(GL_DEPTH_TEST);	//enable Depth Testing

    glutDisplayFunc(display);	//display callback function
    glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occurring)

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMainLoop();		//The main loop of OpenGL

    //capture();

    return 0;
}