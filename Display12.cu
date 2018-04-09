#define _USE_MATH_DEFINES

#include "Display12.h"

#include <iostream>

void displayFunc();
void init();
void reshape(int w,int h);
void idle();

float getMaxPhero(){
    float max = -1;
    for (int i=0; i<MACRO_MAX; i++){
        for (int j=0; j<MACRO_MAX; j++){
            float p = cells[i][j].phero;
            if (max<p){
                max=p;
            }
        }
    }
    return max;
}


void getHeatMapColor(float value, float cl[3]){
    const int NUM_COLORS = 4;
    static float color[NUM_COLORS][3] = { {0,0,1}, {0,1,0}, {1,1,0}, {1,0,0} };
    // A static array of 4 colors:  (blue,   green,  yellow,  red) using {r,g,b} for each.

    int idx1;        // |-- Our desired color will be between these two indexes in "color".
    int idx2;        // |
    float fractBetween = 0;  // Fraction between "idx1" and "idx2" where our value is.

    if(value <= 0){
        idx1 = idx2 = 0;
    }    // accounts for an input <=0
    else if(value >= 1){
        idx1 = idx2 = NUM_COLORS-1;
    }    // accounts for an input >=0
    else{
        value = value * (NUM_COLORS-1);        // Will multiply value by 3.
        idx1  = floor(value);                  // Our desired color will be after this index.
        idx2  = idx1+1;                        // ... and before this index (inclusive).
        fractBetween = value - float(idx1);    // Distance between the two indexes (0-1).
    }

    cl[0] = (color[idx2][0] - color[idx1][0])*fractBetween + color[idx1][0];
    cl[1] = (color[idx2][1] - color[idx1][1])*fractBetween + color[idx1][1];
    cl[2] = (color[idx2][2] - color[idx1][2])*fractBetween + color[idx1][2];
}

void display(int argc, char *argv[]){
    std::cout << "Disp Init" << std::endl;
    glutInitWindowSize(600, 600);
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA | GLUT_DOUBLE);
    glutCreateWindow("Ants");
    glutDisplayFunc(displayFunc);
    glutIdleFunc(idle);
    glutReshapeFunc(reshape);
    init();
    glutMainLoop();
    std::cout << "Disp End" << std::endl;
}

void idle(void)
{
    glutPostRedisplay();
}

void glDrawHexLine(double x,double y, double r){
    glBegin(GL_LINE_LOOP);
    glColor3f(1.0, 1.0, 1.0);
    for (int i=0; i<6; i++){
        glVertex2f(x + r*cos(i*M_PI/3), y + r*sin(i*M_PI/3));
    }
    glEnd();
}

void glDrawHex(double x,double y, double r, float rgb[3]){
    glBegin(GL_TRIANGLE_FAN);
    glColor3f(rgb[0], rgb[1], rgb[2]);
    glVertex2f(x, y);
    for (int i=0; i<=6; i++){
        glVertex2f(x + r*cos(i*M_PI/3), y + r*sin(i*M_PI/3));
    }
    glEnd();
}

void glDrawNest(double x,double y,double r){
    glBegin(GL_TRIANGLE_FAN);
    glColor3f(0.0, 1.0, 0.0);
    glVertex2f(x, y);
    for (int i=0; i<=6; i++){
        glVertex2f(x + r*cos(i*M_PI/3), y + r*sin(i*M_PI/3));
    }
    glEnd();
}

void glDrawFood(double x,double y,double r){
    glBegin(GL_TRIANGLE_FAN);
    glColor3f(1.0, 0.0, 0.0);
    glVertex2f(x, y);
    for (int i=0; i<=6; i++){
        glVertex2f(x + r*cos(i*M_PI/3), y + r*sin(i*M_PI/3));
    }
    glEnd();
}

void drawCells(){
    float maxPhero = getMaxPhero();
    if (maxPhero <= MACRO_EPS ){
        maxPhero = 1.0;
    }
    for(int i=1; i<MACRO_MAX; i++){
        for(int j=1; j<MACRO_MAX; j++){
            double x = cells[i][j].cart.x;
            double y = cells[i][j].cart.y;
            float p = cells[i][j].phero;
            float rgb[3];
            if ( (cells[i][j].status&NEAR_NEST) != NORMAL_CELL ){
                glDrawNest(x,y,1/sqrt(3));
            }
            else if ( (cells[i][j].status&NEAR_FOOD) != NORMAL_CELL ){
                glDrawFood(x,y,1/sqrt(3));
            }
            else{
                getHeatMapColor(p/maxPhero, rgb);
                glDrawHex(x,y,1/sqrt(3),rgb);
            }
        }
    }
}

void drawAnts(){
    for(int k=1; k<MACRO_NMAX; k++){
        int i,j;
        enum AntStatus s = ants[k].status;
        i = ants[k].i;
        j = ants[k].j;
        double x,y;
        x = cells[i][j].cart.x;
        y = cells[i][j].cart.y;

        glPointSize(5);
        glBegin(GL_POINTS);
        switch (s){
            case FORAGE:
                glColor3f(0.0, 1.0, 0.0);
                break;
            case GOHOME:
                glColor3f(1.0, 0.0, 0.0);
                break;
            case EMERGENCY:
                glColor3f(1.0, 1.0, 1.0);
                break;
            default:
                exit(1);
        }
        glVertex2f(x , y);
        glEnd();
    }
}

void displayFunc(void){
    glClear(GL_COLOR_BUFFER_BIT);
    calculation();
    cudaMemcpyFromSymbol(cells,cells_d,MACRO_MAX*MACRO_MAX*sizeof(Cell),0);
    cudaMemcpyFromSymbol(ants,ants_d,MACRO_NMAX*sizeof(Ant),0);
    /* drawCells(); */
    drawAnts();

    glutSwapBuffers();
}

void init(){
    glClearColor(0.0, 0.0, 0.0, 1.0);
}

void reshape(int w, int h){
    glMatrixMode(GL_PROJECTION);
    glViewport(0, 0, w, h);
    glLoadIdentity();
    gluOrtho2D(-(MACRO_CART_X_ZERO+1),(MACRO_MAX-MACRO_CART_X_ZERO)+1,-(MACRO_CART_Y_ZERO+1),(MACRO_MAX-MACRO_CART_Y_ZERO)+1);
}
