#include "MainFrame.h"
#include <iostream>

namespace {

float scale = 1.f;
float aspect = 1.f;

#ifdef __APPLE__
unsigned int SCR_WIDTH = 600;
unsigned int SCR_HEIGHT = 600;
#else
unsigned int SCR_WIDTH = 1000;
unsigned int SCR_HEIGHT = 1000;
#endif

void ScrollCallback(GLFWwindow* window, double xoffset, double yoffset) {
    scale *= std::pow(1.1f, (float)yoffset);
}

void FrameBufferSizeCallback(GLFWwindow* window, int width, int height) {
    SCR_WIDTH = width;
    SCR_HEIGHT = height;
    glViewport(0, 0, SCR_WIDTH, SCR_HEIGHT);        // Set the viewport to cover the new window                  
    aspect = (float)SCR_WIDTH / (float)SCR_HEIGHT;  // Set the aspect ratio of the clipping area to match the viewport
}

}

void MainFrame::LeftMouseMove(float start_x, float start_y, float curr_x, float curr_y) {
    if (modeling_state_ == OBJ_ROTATION) {
        // ---------------------------------- Object Rotation ---------------------------------------
        // TODO: Add your code here.
        // Find the correct 4x4 transform matrix "transform_mat" to rotate the object about its center.

        glm::mat4x4 transform_mat(1.f);

        glm::vec2 V(curr_x - start_x, curr_y - start_y);
        glm::vec2 A(-V.y,V.x);

        std::tuple<glm::vec3, glm::vec3> screen_world_start_ray = Screen2WorldRay(start_x, start_y);

        std::tuple<int, glm::vec3> intersection_tuple_start = mesh_.FaceIntersection(std::get<0>(screen_world_start_ray), std::get<1>(screen_world_start_ray));

        glm::vec3 Sstart_world = Screen2World(start_x, start_y, std::get<1>(intersection_tuple_start).z);
        glm::vec3 Sstart_A_world = Screen2World(start_x + A.x, start_y + A.y, std::get<1>(intersection_tuple_start).z);

        glm::vec3 ra = glm::normalize(Sstart_A_world-Sstart_world);
        float angle = 0.01f * std::sqrt(std::pow(A.x,2)+std::pow(A.y,2));
        



        glm::mat4x4 R;
        R = glm::rotate(glm::mat4x4(1.f), angle, ra);

        glm::mat4x4 T({1.f,0.0f,0.0f,0.0f, 0.0f,1.f,0.0f,0.0f, 0.0f,0.0f,1.f,0.0f, -mesh_.center_.x,-mesh_.center_.y,-mesh_.center_.z,1.f});
        glm::mat4x4 T_inv({1.f,0.0f,0.0f,0.0f, 0.0f,1.f,0.0f,0.0f, 0.0f,0.0f,1.f,0.0f, mesh_.center_.x,mesh_.center_.y,mesh_.center_.z,1.f});


        transform_mat = T_inv * R * T;
        // ...
        
        transform_mat = glm::transpose(transform_mat);


        mesh_.ApplyTransform(transform_mat);
    }
    else if (modeling_state_ == OBJ_TRANSLATION) {
        // ---------------------------------- Object Translation ------------------------------------
        // TODO: Add your code here.
        // Find the correct 4x4 transform matrix "trans_mat" to translate the object along the view plane.

        glm::mat4x4 transform_mat(1.f);

        

        std::tuple<glm::vec3, glm::vec3> screen_world_start_ray = Screen2WorldRay(start_x, start_y);
        // std::tuple<glm::vec3, glm::vec3> screen_world_curr_ray = Screen2WorldRay(curr_x, curr_y);

        std::tuple<int, glm::vec3> intersection_tuple_start = mesh_.FaceIntersection(std::get<0>(screen_world_start_ray), std::get<1>(screen_world_start_ray));

        glm::vec3 Pstart = std::get<1>(intersection_tuple_start);
        
        glm::vec3 Pcurr = Screen2World(curr_x, 
                                        curr_y, 
                                        std::get<1>(screen_world_start_ray).z-std::get<0>(screen_world_start_ray).z); 

        // calculate Pstart again because some bug on function Screen2World()
        Pstart = Screen2World(start_x, 
                            start_y, 
                            std::get<1>(screen_world_start_ray).z-std::get<0>(screen_world_start_ray).z); 

        
        transform_mat = glm::translate(transform_mat, Pcurr - Pstart);

        




        mesh_.ApplyTransform(transform_mat);
        
    }
    else if (modeling_state_ == OBJ_EXTRUDE) {
        // ---------------------------------- Face Extrusion ------------------------------------
        // TODO: Add your code here.
        // Find the correct 4x4 transform matrix "trans_mat" to translate the face vertices along the face normal.
 
        glm::mat4x4 transform_mat(1.f);
        int face_index = -1;

        // glm::vec3 cam_pos = camera_.view_mat_inv_ * glm::vec4(0.f, 0.f, 0.f, 1.f);

        // get the ray R0 which goes through Pstart
        std::tuple<glm::vec3, glm::vec3> screen_world_start_ray = Screen2WorldRay(start_x, start_y);

        // get the ray R1 which goes through Pcurr
        std::tuple<glm::vec3, glm::vec3> screen_world_curr_ray = Screen2WorldRay(curr_x, curr_y);
        glm::vec3 R0 = std::get<1>(screen_world_start_ray);

        // get the <int>face_index and <vec3>hit_point
        std::tuple<int, glm::vec3> intersection_tuple_start = mesh_.FaceIntersection(std::get<0>(screen_world_start_ray), R0);
        glm::vec3 Pstart =  std::get<1>(intersection_tuple_start);
        face_index = std::get<0>(intersection_tuple_start);

        // normal vector of M,V,Q
        glm::vec3 M_normal = glm::cross( mesh_.vertices_[mesh_.faces_[face_index][0]] - Pstart, mesh_.vertices_[mesh_.faces_[face_index][1]] - Pstart);
        glm::vec3 V = glm::cross(M_normal,std::get<1>(screen_world_curr_ray));
        glm::vec3 Q_normal = glm::cross(V,std::get<1>(screen_world_curr_ray));

        // projection: project R0 on Q
        glm::vec3 R0_Q = (glm::dot(R0,Q_normal)/(Q_normal/glm::normalize(Q_normal)));

        // projection: project R0_Q on M
        glm::vec3 t = (glm::dot(R0_Q,M_normal)/(M_normal/glm::normalize(M_normal)));

        // Pcurr
        glm::vec3 Pcurr = Pstart + t * M_normal;
        
        transform_mat = glm::translate(transform_mat, Pcurr - Pstart);
        
        // ...

        mesh_.ApplyFaceTransform(face_index, transform_mat);
    }
}

void MainFrame::VisualizeWorldSpace() {
    // ---------------------------------- World Space Visualization ------------------------------------
    // TODO: Add your code here to visualize the world space.


    glBegin(GL_LINES);
    glColor3f(1.f,0.0f,0.0f);
    glVertex3f(0.0f,0.0f,0.0f);
    glColor3f(1.f,0.0f,0.0f);
    glVertex3f(2.0f,0.0f,0.0f);
    glColor3f(0.f,1.f,0.f);
    glVertex3f(0.0f,0.0f,0.0f);
    glColor3f(0.f,1.f,0.f);
    glVertex3f(0.0f,2.0f,0.0f);
    glColor3f(0.f,0.f,1.f);
    glVertex3f(0.0f,0.0f,0.0f);
    glColor3f(0.f,0.f,1.f);
    glVertex3f(0.0f,0.0f,2.0f);
    glEnd();


    


    // ...
}


// -------------------------------------------------------------------------------------
// -------------------------- No need to change ----------------------------------------
// -------------------------------------------------------------------------------------
void MainFrame::MainLoop() {
    // glfw: initialize and configure
    glfwInit();
    glfwWindowHint(GLFW_SAMPLES, 4);

    // glfw window creation, set viewport with width=1000 and height=1000
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "3DModeling", NULL, NULL);
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, FrameBufferSizeCallback);
    glfwSetScrollCallback(window, ScrollCallback);
    // glad: load all OpenGL function pointers
    gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);

    const float alpha = 0.3f;
    const float beta = 0.1f;

    const float r = 5.f;
    camera_.LookAt(r * glm::vec3(std::cos(alpha) * std::cos(beta), std::cos(alpha) * std::sin(beta), std::sin(alpha)),
        glm::vec3(0.f, 0.f, 0.f),
        glm::vec3(0.f, 0.f, 1.f));

    glEnable(GL_DEPTH_TEST);

    // render loop
    while (!glfwWindowShouldClose(window)) {
        ProcessInput(window);

        // glEnable(GL_DEPTH_TEST);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();

        // Apply camera projection;
        camera_.Perspective(90.f, aspect, .5f, 10.f);
        camera_.UpdateScale(scale);
        scale = 1.f;
        camera_.ApplyProjection();

        glClearColor(0.0, 0.0, 0.0, 1.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);        // Clear the display

        DrawScene();

        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        glfwPollEvents();
        glfwSwapBuffers(window);
    }

    // glfw: terminate, clearing addl previously allocated GLFW resources.
    glfwTerminate();
}

void MainFrame::ProcessInput(GLFWwindow* window) {
    // Key events
    if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS) {
        modeling_state_ = OBJ_ROTATION;
    }
    if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS) {
        modeling_state_ = OBJ_TRANSLATION;
    }
    if (glfwGetKey(window, GLFW_KEY_3) == GLFW_PRESS) {
        modeling_state_ = OBJ_SUBDIVIDE;
    }
    if (glfwGetKey(window, GLFW_KEY_4) == GLFW_PRESS) {
        modeling_state_ = OBJ_EXTRUDE;
    }

    int current_l_mouse_state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);

    // Handle left mouse
    if (current_l_mouse_state == GLFW_PRESS) {
        double xposd, yposd;
        float xpos, ypos;
        glfwGetCursorPos(window, &xposd, &yposd);
        xpos = float(xposd);
        ypos = float(SCR_HEIGHT - yposd);
        if (l_mouse_state_ == GLFW_RELEASE) {
            LeftMouseClick(xpos, ypos);
            l_click_cursor_x_ = xpos;
            l_click_cursor_y_ = ypos;
        }
        if (l_mouse_state_ == GLFW_PRESS &&
            (std::abs(xpos - last_cursor_x_) > 2.f || std::abs(ypos - last_cursor_y_) > 2.f)) {
            LeftMouseMove(l_click_cursor_x_, l_click_cursor_y_, xpos, ypos);
        }
        last_cursor_x_ = float(xpos);
        last_cursor_y_ = float(ypos);
    }
    if (current_l_mouse_state == GLFW_RELEASE) {
        if (l_mouse_state_ == GLFW_PRESS) {
            LeftMouseRelease();
        }
    }
    l_mouse_state_ = current_l_mouse_state;

    // Handle right mouse
    int current_r_mouse_state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT);
    if (current_r_mouse_state == GLFW_PRESS) {
        double xposd, yposd;
        float xpos, ypos;
        glfwGetCursorPos(window, &xposd, &yposd);
        xpos = float(xposd);
        ypos = float(SCR_HEIGHT - yposd);
        if (r_mouse_state_ == GLFW_RELEASE) {
            RightMouseClick(xpos, ypos);
        }
        if (r_mouse_state_ == GLFW_PRESS &&
            (std::abs(xpos - last_cursor_x_) > 2.f || std::abs(ypos - last_cursor_y_) > 2.f)) {
            RightMouseMove(last_cursor_x_, last_cursor_y_, xpos, ypos);
        }
        last_cursor_x_ = float(xpos);
        last_cursor_y_ = float(ypos);
    }
    if (current_r_mouse_state == GLFW_RELEASE) {
        if (r_mouse_state_ == GLFW_PRESS) {
            RightMouseRelease();
        }
    }
    r_mouse_state_ = current_r_mouse_state;
}

void MainFrame::LeftMouseClick(float x, float y) {
    if (modeling_state_ == OBJ_SUBDIVIDE) {
        glm::vec3 p_world = Screen2World(x, y);
        glm::vec3 cam_pos = camera_.view_mat_inv_ * glm::vec4(0.f, 0.f, 0.f, 1.f);
        mesh_.SubdivideFace(cam_pos, glm::normalize(p_world - cam_pos));
    }
    else if (modeling_state_ == OBJ_EXTRUDE) {
        glm::vec3 p_world = Screen2World(x, y);
        glm::vec3 cam_pos = camera_.view_mat_inv_ * glm::vec4(0.f, 0.f, 0.f, 1.f);
        mesh_.GenExtrudeFace(cam_pos, glm::normalize(p_world - cam_pos));
    }
}

void MainFrame::LeftMouseRelease() {
    mesh_.CommitTransform();
}

void MainFrame::RightMouseClick(float x, float y) {
    return;
}

void MainFrame::RightMouseMove(float start_x, float start_y, float curr_x, float curr_y) {
    glm::vec2 s_start(start_x, start_y);
    glm::vec2 s_cur(curr_x, curr_y);
    glm::vec2 V = s_cur - s_start;
    glm::vec2 A = glm::vec2(-V.y, V.x);
    glm::vec3 rot_axis = glm::normalize(Screen2World(A + s_start) - Screen2World(s_start));
    glm::mat4x4 rot_mat = glm::rotate(glm::mat4x4(1.f), 0.007f * glm::length(A), rot_axis);
    camera_.ApplyTransform(rot_mat);
}

void MainFrame::RightMouseRelease() {
    return;
}

glm::vec3 MainFrame::Camera2World(const glm::vec3& x, float w) {
    return glm::vec3(camera_.view_mat_inv_ * glm::vec4(x, w));
}

glm::vec3 MainFrame::World2Camera(const glm::vec3& x, float w) {
    return glm::vec3(camera_.view_mat_ * glm::vec4(x, w));
}

glm::vec3 MainFrame::Screen2World(const glm::vec2& v, float depth) {
    float x = v.x / SCR_WIDTH  * 2.f - 1.f;
    float y = v.y / SCR_HEIGHT * 2.f - 1.f;
    float focal = std::tan(camera_.fov_ * .5f / 180.f * glm::pi<float>());
    glm::vec4 v_camera(x * focal * aspect, y * focal, -1.f, 1.f);
    v_camera = v_camera * depth;
    glm::vec4 v_world = camera_.view_mat_inv_ * v_camera;
    return glm::vec3(v_world);
}

glm::vec3 MainFrame::Screen2World(float scr_x, float scr_y, float camera_z) {
    float x = scr_x / SCR_WIDTH * 2.f - 1.f;
    float y = scr_y / SCR_HEIGHT * 2.f - 1.f;
    float focal = std::tan(camera_.fov_ * .5f / 180.f * glm::pi<float>());
    glm::vec4 v_camera(x * focal * aspect, y * focal, -1.f, 1.f);
    v_camera = v_camera * -camera_z;
    glm::vec4 v_world = camera_.view_mat_inv_ * v_camera;
    return glm::vec3(v_world);
}

std::tuple<glm::vec3, glm::vec3> MainFrame::Screen2WorldRay(float scr_x, float scr_y) {
    float x = scr_x / SCR_WIDTH * 2.f - 1.f;
    float y = scr_y / SCR_HEIGHT * 2.f - 1.f;
    float focal = std::tan(camera_.fov_ * .5f / 180.f * glm::pi<float>());
    glm::vec3 o = camera_.view_mat_inv_ * glm::vec4(0.f, 0.f, 0.f, 1.f);
    glm::vec4 v_camera(x * focal * aspect, y * focal, -1.f, 0.f);
    glm::vec3 v = camera_.view_mat_inv_ * v_camera;
    return std::make_tuple(o, v);
}

void MainFrame::DrawScene() {
    // Draw mesh
    mesh_.Draw();

    VisualizeWorldSpace();
}