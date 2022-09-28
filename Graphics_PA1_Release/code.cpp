#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "code.h"
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <vector>


//////////////////////////////////////////////////////
// The Global Variables To Be Used

struct Triangle {
    double vertices[3][2];
    double matrix[3][3];
    int color_index;
};

struct Matrix {
    double matrix[3][3];
};

// list of triangles, each element is of type Triangle
// to access the kth triangle, just use triangles[k]
// to get the number of triangles in the list, use triangles.size()
std::vector<Triangle> triangles;

// a triangle object for temporary storage of points
Triangle triangle_to_draw;
// count the number of points specified in this triangle
int point_count = 0;

// depth of recursion for IFS
// inital value is 8
// change the initial value here
int recursion_depth = 8;


// color array for triangles
// size is 11. So color_index should range from 0 to 10 for triangles.
double color_array[][3] = {
    {0.9,0,0}, //red
    {0,0.5,0.4},
    {0.1,0.2,0.46},
    {0.9,0.9,0},
    {0,1.0,0},
    {0,1.0,1.0},
    {0,0,1.0},
    {1.0,0,1.0},
    {0.9,0.6,0},
    {0.9,1.0,0.6},
    {0.2,0.2,0.2}
};

std::vector<Matrix> affine_matrices;

void VerticesToMatrix(double vertices[][2], double matrix[][3]){
    for ( int i = 0; i < 3; i++ ){
        for ( int j = 0; j < 2; j++){
            matrix[j][i] =  vertices[i][j];
        }
    }
    for ( int i = 0; i < 3; i++ ){
        matrix[2][i] = 1.0;
    }
}

void MatrixMultiplication(double matrix_1[][3], double matrix_2[][3], double result[][3]){
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            result[i][j] = 0;
 
            for (int k = 0; k < 3; k++) {
                result[i][j] += matrix_1[i][k] * matrix_2[k][j];
            }
        }

    }
}

void PrintMatrix(double matrix[3][3]){
    for (int i=0; i<3;i++){
        for (int j=0; j<3;j++){
            std::cout << matrix[i][j] << " " << std::ends;
        }
        std::cout << "\n" << std::ends;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////
// This function is called to clear triangles.

void ClearTriangles() {
    triangles.clear();
    point_count = 0;
}

/////////////////////////////////////////////////////////////////////////////////////////
// This function is called to draw triangles in the Triangles window.

void DrawTriangles() {
    // uncomment this line if you would like to unfill triangles
    // glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);

    // sample code to draw vertices of triangle
    // glColor3d(1.0, 1.0, 1.0);
    // glPointSize(4);
    // glBegin(GL_POINTS);
    // for (int i = 0; i < point_count; i++) {
    //     glVertex2d(triangle_to_draw.vertices[i][0], triangle_to_draw.vertices[i][1]);
    // }
    // glEnd();

    // TODO: Add code to draw triangles here. Use triangles.size() to get number of triangles.
    // Use triangles[i] to get the ith triangle.
    // Remember to set current gl color from the color_array

    glBegin(GL_TRIANGLES);
    for (int i = 0; i < point_count/3; i++) {
        glColor3d( color_array[i][0], color_array[i][1], color_array[i][2] );
        for (int j = 0; j < 3; j++){
            glVertex2d(triangles[i].vertices[j][0], triangles[i].vertices[j][1]);
        }
            
    }
    glEnd();
    
}


/////////////////////////////////////////////////////////////////////////////////////////
// This function draws the factal in a recursion manner. The depth of recursion is k.
// When k=0, the function draw the most basic shape, that is, the original triangle, and
// stops recursion.

void RecursiveFractal(int k) {

    // TODO: Add code to implement the IFS method to draw fractal.
    // You can follow the pseudo code in the slides.
    // Use triangles.size() to get number of triangles.
    // Use triangles[i] to get the ith triangle.
    // The fields of struct Triangle include:
    //	 double vertices[3][2];
    //	 double matrix[3][3];

    if ( k > 0 ){
        for (int i = 0; i < triangles.size(); i++){
            glPushMatrix();
            glTranslated(-0.25, -0.25, 0.0); 
            // glScaled(0.5, 0.5, 0.0);
            glMultMatrixd((const GLdouble *)affine_matrices[i].matrix);
            RecursiveFractal(k-1);
            // for (int j = k-1; j <= i; j++){
            //     DrawTriangles();
            // }

            glPopMatrix();

            glPushMatrix();
            glTranslated(0.25, -0.25, 0.0); 
            // glScaled(0.5, 0.5, 0.0);
            glMultMatrixd((const GLdouble *)affine_matrices[i].matrix);
            RecursiveFractal(k-1);
            // for (int j = k-1; j <= i; j++){
            //     DrawTriangles();
            // }

            glPopMatrix();

            glPushMatrix();
            glTranslated(0, 0.25, 0.0); 
            // glScaled(0.5, 0.5, 0.0);
            glMultMatrixd((const GLdouble *)affine_matrices[i].matrix);
            RecursiveFractal(k-1);
            // for (int j = k-1; j <= i; j++){
            //     DrawTriangles();
            // }
            glPopMatrix();
        }
    }else{
        glBegin(GL_TRIANGLES);
        glVertex2d(triangles[0].vertices[0][0], triangles[0].vertices[0][1]);
        glVertex2d(triangles[0].vertices[1][0], triangles[0].vertices[1][1]);
        glVertex2d(triangles[0].vertices[2][0], triangles[0].vertices[2][1]);
        glEnd();

    }

}

/////////////////////////////////////////////////////////////////////////////////////////
// This function invokes RecursiveFractal()

void ConstructiveFractals() {

    if (triangles.size() < 2)
        return;


    glColor3f(1.0, 1.0, 0.0);
    RecursiveFractal(recursion_depth);
}

/////////////////////////////////////////////////////////////////////////////////////////
// This function is called to handle mouse left click events.
// m_x and m_y are the x and y coordinates of the clicked point in OpenGL coordinate system.

void MouseInteraction(double m_x, double m_y) {
    point_count += 1;

    std::cout << "x = " << m_x << "; y = " << m_y<< std::endl;
    // TODO: Store the point temporarily into the variable triangle_to_draw.
    triangle_to_draw.vertices[ (point_count-1)%3 ][0] = m_x;
    triangle_to_draw.vertices[ (point_count-1)%3 ][1] = m_y;
    // When 3 points are specified, we get a new triangle.
    if ( point_count % 3 == 0 ){

        if ( point_count >= 6)
            AffineMatricesCalculation(
                triangles[0].vertices, 
                triangle_to_draw.vertices, 
                triangle_to_draw.matrix
                );

        triangles.push_back(triangle_to_draw);
    }
    // Compute the matrix for affine transformation from the first triangle to this one
    //	 by invoking AffineMatricesCalculation().
    // Store both the points and matrix of the triangle into a new element of the list 'triangles'
    // Use push_back function from std::Vector to add new triangles to the list

}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// This function is a tool function that computes matrix for the affine transformation from one original
//   triangle to another triangle.
// Input:
//		 v_original[][2]	the pointer to an array containing data of original triangle
//		 v_transformed[][2]	the pointer to an array containing data of triangle obtained by transforming original triangle
// Output:
//		 matrix[][3]		a pointer to 3x3 matrix that is the affine transformation that changes original triangle to the later one.

void AffineMatricesCalculation(double v_original[][2], double v_transformed[][2], double matrix[][3]) {

    // TODO: Compute the affine transformation matrix that transforms triangle specified in v_original to the one specified in v_transformed.
    //		 Base the computation on the formula M = T'T^(-1), where T' is the 3x3 matrix with each column the homogeneous coordinates of transformed triangle's vertex
    //		 and T is 3x3 matrix organized in a similar manner but stores data of the original triangle.
    //		 If you do not want to calculate the inverse of T yourself, we provide a tool function InverseMatrix(). This function could compute the inverse of T.
    


    Matrix result_matrix;
    // temporarily store the original matrix, to be inversed
    Matrix transformed_m;
    Matrix inversed_m;
    // vertices to matrix -------------------------------------------------------
    VerticesToMatrix(triangles[0].vertices, triangles[0].matrix);
    VerticesToMatrix(triangle_to_draw.vertices, transformed_m.matrix);
    // --------------------------------------------------------------------------
    std::cout << "triangles[0].matrix" << std::endl;
    PrintMatrix(triangles[0].matrix);

    // store matrix to the correspon triangle in the "Triangles" list -----------
    for (int i = 0; i<3; i++){
        for (int j = 0; j<3; j++){

                matrix[i][j] = transformed_m.matrix[i][j];

        }
    }
    // --------------------------------------------------------------------------
    InverseMatrix(triangles[0].matrix, inversed_m.matrix);
    MatrixMultiplication(transformed_m.matrix, triangles[0].matrix, result_matrix.matrix);

    std::cout << "inversed_m.matrix" << std::endl;
    PrintMatrix(inversed_m.matrix);

    std::cout << "result_matrix.matrix" << std::endl;
    PrintMatrix(result_matrix.matrix);

    // no original matrix, only all the transformation matrices
    affine_matrices.push_back(result_matrix);

}


////////////////////////////////////////////////////////////////////////////////////////////////////////
// A routine to calculate inverse matrix of a 3x3 matrix which has all its values in the third row being 1.
//	original_m: 3x3 matrix with original_m[2][0]=original_m[2][1]=original_m[2][2]=1.
//  inverse_m:  3x3 matrix, the inverse of original_m.
//
void InverseMatrix(double original_m[][3], double inverse_m[][3]) {
    double determinant;
    determinant = original_m[0][0] * (original_m[1][1] - original_m[1][2]) - original_m[0][1] * (original_m[1][0] - original_m[1][2]) +
                  original_m[0][2] * (original_m[1][0] - original_m[1][1]);

    inverse_m[0][0] = (original_m[1][1] - original_m[1][2]) / determinant;
    inverse_m[1][0] = (original_m[1][2] - original_m[1][0]) / determinant;
    inverse_m[2][0] = (original_m[1][0] - original_m[1][1]) / determinant;

    inverse_m[0][1] = (original_m[0][2] - original_m[0][1]) / determinant;
    inverse_m[1][1] = (original_m[0][0] - original_m[0][2]) / determinant;
    inverse_m[2][1] = (original_m[0][1] - original_m[0][0]) / determinant;

    inverse_m[0][2] = (original_m[0][1] * original_m[1][2] - original_m[0][2] * original_m[1][1]) / determinant;
    inverse_m[1][2] = (original_m[0][2] * original_m[1][0] - original_m[0][0] * original_m[1][2]) / determinant;
    inverse_m[2][2] = (original_m[0][0] * original_m[1][1] - original_m[0][1] * original_m[1][0]) / determinant;
}