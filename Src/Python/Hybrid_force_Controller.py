import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
from scipy.integrate import solve_ivp
import cvxpy as cp
from mpl_toolkits.mplot3d.art3d import Poly3DCollection



#Initializing variables

i_goal = 0
time_start = None
time_prev = None
goal_threshold = 0.5
force_switch = 0
delta_uf_prev = np.zeros(3)

Projection_matrix = np.zeros((3,3))

F_Measured = np.zeros(3)

error = np.zeros(3)
delta_ux_prev = np.zeros(3)

c = 0


F_sensor_reading = np.zeros(3)





i_goal = 0



stage = 0 # 0: when SBA is trying to reach the goal_point
          # 1: goal_point is reached and force_control starts and goes till it maintains the desired force for desired time
          # 2: desired force for desired duration is complete and SBA starts to deflate till no force is exerted

on_plane = False #False if SBA's tip hasn't reached the plane
                 #True  if tip is on the plane


DT = 0.1

#Initializing the spatial_error integral
Ex_current  = np.zeros(3)
Ef_current  = np.zeros(3)



class Kinematics:
    def __init__(self):

        self.n    = 3                 #Number of segments for the SBA
        self.d    = 4.7               #Geometric parameter
        self.gain = 0.5                 #Gain for the U_ell
        self.r    = 5
        self.a    = 2
        self.h    = 0.0381
        self.E    = 7.35
        self.v    = 0.5
        self.K    = np.diag([15,15,5,15,15,5,15,15,5])
        self.D    = np.diag([5,5,5,5,5,5,5,5,5])
        self.I    = np.diag((self.K.diagonal()**2)/(4*self.D.diagonal()))
        self.A    = 25.844e-6
        self.Kf   = np.diag([0.5,0.5,0.5])
        self.dl   = np.diag([5,5,5])
        self.kl   = np.diag([.1,.1,.1])
        self.K_f  = 1
        self.I_f  = 2
        self.theta_plane    = -np.pi/4  #orientation of plane w.r.t global frame
        self.P_plane_origin = np.array([-9.67,0,27.02]) #coordinates of plane center w.r.t global frame
        self.sensor_var = 0.009  #0.009 N , assuming axes are not correlated and each have same variance

    def signM(self, x):
        sign_x = x/np.abs(x)
        return sign_x


    def pressureToVolume(self, pressures):
        '''In: pressues in kpa
        out: volume as a fucntion of pressures in mL
        '''
        # print(pressures)
        # print(pressures[0])

        v1 = self.signM(pressures[0])*0.161*np.abs(pressures[0])**(0.333)
        v2 = self.signM(pressures[1])*0.161*np.abs(pressures[1])**(0.333)
        v3 = self.signM(pressures[2])*0.161*np.abs(pressures[2])**(0.333)

        return np.array([v1, v2, v3])


    def pressureToLength(self, P, linspaceParam):
        # Use input pressure to determine chamber lengths based on timoshenko's model
        # corrected model for fresh SBA, scale = 0.84
        scale = 0.84
        len_c = scale*2*linspaceParam*0.622*self.a*( P*self.a/(self.E*self.h) )**(1/3)

        return len_c

    def volume_to_length(self, V):
        l1, l2, l3 = 45.67*V[0],  45.67*V[1],  45.67*V[2]

        return np.array([l1, l2, l3])

    def lengths_to_volumes(self, ells, base_height):
        volumes = 100/4667 * (ells - base_height)
        return volumes*1000

    def Force_sensor(self, x_no_load, x_model):

        #Measures the force from  the environment on the robot

        F = self.Kf@(x_no_load - x_model)
        return F

    def Rz(self, phi):
        T = np.array([[np.cos(phi), -np.sin(phi), 0],
                      [np.sin(phi),  np.cos(phi), 0],
                      [          0,            0, 1]])
        return T

    def Ry(self, theta):
        T = np.array([[ np.cos(theta), 0, np.sin(theta)],
                      [             0, 1,             0],
                      [-np.sin(theta), 0, np.cos(theta)]])
        return T

    def Rx(self, theta):
        T = np.array([[1,             0,              0],
                      [0, np.cos(theta), -np.sin(theta)],
                      [0, np.sin(theta),  np.cos(theta)]])
        return T


    def Link_function(self, qi): #Method to calculate the length of each prismatic links of each segment
        theta_i = qi[1,0]        #based on current configuration variables
        S_i     = qi[2,0]        #also evaluating the variation of these link lengths w.r.t change in theta, and arc length
                                 #The link length does not depend upon the out of the plane angle (phi)
        if theta_i < 1e-3:
            li = S_i
            dl_theta = 0
            dl_S = 1
        else:
            li = 2*S_i*np.sin(theta_i/2)/theta_i
            dl_theta = 2*S_i*(np.cos(theta_i/2)/(2*theta_i) - np.sin(theta_i/2)/(theta_i**2))
            dl_S = 2*np.sin(theta_i/2)/theta_i

        return li, dl_theta, dl_S

    def single_map(self, qi):
        #Method helps to evaluate the Xi vector along with the Jacobian_xi_w.r.t q for each segment
        phi   = qi[0,0]
        theta = qi[1,0]

        li, dl_theta, dl_S = self.Link_function(qi)

        xi = np.array([[phi], [theta/2], [li], [theta/2], [-phi]])

        Jmi = np.array([[1,0,0], [0, 0.5, 0], [0,dl_theta,dl_S], [0,0.5,0],[-1,0,0]])

        return xi, Jmi

    def Build_XI(self, q):
        #Method helps to build the final Xi which is a (3*n,1) vector and J_xi_q which is a (5*n, 3*n)
        xi = np.zeros((5*self.n,1))
        Jm = np.zeros((5*self.n, 3*self.n))

        for i in range(self.n):
            rows = slice(i*5, i*5 + 5)
            qi = q[3*i : 3*i + 3].reshape(3,1)
            xi_block, Jmi = self.single_map(qi)

            xi[rows] = xi_block
            Jm[rows, 3*i : 3*i + 3] = Jmi
        return xi, Jm

    def Augmented_DH(self, q):
        #Function to derive the DH Table based on our revolute and prismatic joints
        DH_aug = np.zeros((5*self.n, 4))
        Xi, _ = self.Build_XI(q)

        for i in range(self.n):
            # if i == (self.n - 1):
            #     base = i*5
            #     xi = Xi[base: base + 5, 0]
            #     DH_aug[base + 0, :] = [xi[0],    0,    0, -np.pi/2]
            #     DH_aug[base + 1, :] = [xi[1],    0,    0,  np.pi/2]
            #     DH_aug[base + 2, :] = [0,       xi[2], 0, -np.pi/2]
            #     DH_aug[base + 3, :] = [xi[3],    0,    0,  np.pi/2]
            #     DH_aug[base + 4, :] = [xi[4],    r,    0,        0]
            # else:
            base = i*5
            xi = Xi[base: base + 5, 0]
            DH_aug[base + 0, :] = [xi[0],    0,    0, -np.pi/2]
            DH_aug[base + 1, :] = [xi[1],    0,    0,  np.pi/2]
            DH_aug[base + 2, :] = [0,       xi[2], 0, -np.pi/2]
            DH_aug[base + 3, :] = [xi[3],    0,    0,  np.pi/2]
            DH_aug[base + 4, :] = [xi[4],    0,    0,        0]

        return DH_aug

    def DH_transform(self, theta, d, a, alpha):
        #Basic 4x4 Tranformation matrix
        T = np.array([[np.cos(theta), -np.sin(theta)*np.cos(alpha),  np.sin(theta)*np.sin(alpha), a*np.cos(theta)],
                      [np.sin(theta),  np.cos(theta)*np.cos(alpha), -np.cos(theta)*np.sin(alpha), a*np.sin(theta)],
                      [0,                            np.sin(alpha),                np.cos(alpha),               d],
                      [0,                                     0,                               0,               1]])

        return T

    def Compute_Augmented_Jacobian(self, dh_aug):

        # This function allows us to evaluate the Jacobian for the RRPRR model. The jacobian size is
        # 3 x 5*n, (here n is number of segments and 5 basically is the number of joints per segment
        # in our case we have 5 joints per segment)

        N = 5*self.n
        T = np.eye(4)
        p = np.zeros((3, N+1))
        z = np.tile(np.array([[0.],[0.],[1.]]), (1 , N+1))

        new_dh_aug = dh_aug
        new_dh_aug[-1,1] = self.r

        for k in range(N):
            theta, d, a, alpha = new_dh_aug[k]
            T     = T@self.DH_transform(theta,d,a,alpha)
            p[:, k+1] = T[0:3, 3]
            z[:, k+1] = T[0:3, 2]

        p_tip = p[:, -1]
        Jxi   = np.zeros((3, N))

        for k in range(N):
            idx = k%5
            if idx == 2:
                Jxi[:,k] = z[0:3, k]
            else:
                Jxi[:, k] = np.cross(z[0:3,k], (p_tip - p[:,k]))

        return Jxi

    def ComputeJacobian(self, q):

        # Evaluating the final jacobian for the soft robot J(q). Its size is 3 x 3*n

        dh_aug = self.Augmented_DH(q)
        Jxi    = self.Compute_Augmented_Jacobian(dh_aug)
        _, Jm  = self.Build_XI(q)

        J_q    = Jxi @ Jm

        return J_q

    def q_diff(self, q1, q2):
        def angular_difference(phi1,phi2):
            diff = (phi1 - phi2 + math.pi) %(2*math.pi) - math.pi
            return diff

        diff_q = np.zeros(3*self.n)
        for i in range(self.n):
            idx = 3*i
            diff_q[idx] = angular_difference(q1[idx], q2[idx])
            diff_q[idx + 1] = q1[idx + 1] - q2[idx + 1]
            diff_q[idx + 2] = q1[idx + 2] - q2[idx + 2]

        return diff_q

    def Model(self, t, q, f_ext_fun, q_ref, Ei_prev):

        # Our basic dynamic model K(q - q_ref) + Dq_dot = J^T.F_ext

        K = self.K
        D = self.D
        I = self.I
        J  = self.ComputeJacobian(q)
        f  = f_ext_fun(t)
        diff_q_q_ref = self.q_diff(q, q_ref)

        Residual = J.T.dot(f) - K.dot(diff_q_q_ref) - I.dot(Ei_prev)
        dq = np.linalg.solve(D, Residual)
        return dq

    def new_Model(self, t, q, f_ext_fun, q_no_load, q_no_load_dot):

        # Our basic dynamic model Kq + D(q_dot - q_ref_dot) = J^T.F_ext

        K = self.K
        D = self.D
        J  = self.ComputeJacobian(q)
        f  = f_ext_fun(t)
        diff_q_no_load_q = self.q_diff(q_no_load, q)
        Residual = J.T.dot(f) + K.dot(diff_q_no_load_q)
        dq = q_no_load_dot + np.linalg.solve(D, Residual)
        return dq

    def q_dynamics(self, q0, Ei_prev, f_ext_fun, t0, dt):

        # Function to evaluate the configuration variable (q) at t = t + dt under the action of
        # external force acting at the tip

        q0 = np.asarray(q0).reshape(-1)
        t_span = (t0, t0 + dt)

        sol = solve_ivp(fun = lambda t, q: self.Model(t, q, f_ext_fun, q0, Ei_prev),
                        t_span = t_span,
                        y0 = q0,
                        method= 'RK45',
                        t_eval = [t0 + dt])

        q_next = sol.y[:,-1]

        return q_next.reshape(-1,1)

    def q_dynamics_new(self, q_init, q_no_load, q_no_load_dot, f_ext_fun, t0, dt):

        # Function to evaluate the configuration variable (q) at t = t + dt under the action of
        # external force acting at the tip

        q_init = np.asarray(q_init).reshape(-1)
        t_span = (t0, t0 + dt)

        sol = solve_ivp(fun = lambda t, q: self.new_Model(t, q, f_ext_fun, q_no_load, q_no_load_dot),
                        t_span = t_span,
                        y0 = q_init,
                        method= 'RK45',
                        t_eval = [t0 + dt])

        q_next = sol.y[:,-1]

        return q_next.reshape(-1,1)

    def Model_dynamics(self,t, q, f_sensor, q_no_load, q_no_load_dot):

        #ODE:  q̇ = q̇_no_load + D⁻¹ [ Jᵀ F + K (q_no_load - q) ].
        #*t* is unused but required by solve_ivp.
        #load is in grams

        K = self.K
        D = self.D
        J = self.ComputeJacobian(q)
        diff_q_no_load_q = self.q_diff(q_no_load, q)
        dq = q_no_load_dot + np.linalg.solve(D, J.T.dot(f_sensor) + K.dot(diff_q_no_load_q))

        return dq

    def q_model_dynamics(self, q_init, q_no_load, q_no_load_dot, f_sensor, t0, dt):
        q_init = np.asarray(q_init).reshape(-1)
        t_span = (t0, t0+dt)

        sol = solve_ivp(fun = lambda t, q: self.Model_dynamics(t,q,f_sensor, q_no_load, q_no_load_dot),
                        t_span= t_span,
                        y0 = q_init,
                        method = 'RK45',
                        t_eval = [t0 + dt])

        q_next_timestep = sol.y[:,-1]

        return q_next_timestep

    def q_no_load(self, ell):

        # Function to evaluate the q under no load from the given Ells

        q_0 = np.zeros(3*self.n)


        ell_1, ell_2, ell_3 = ell

        S = (ell_1 + ell_2 + ell_3)/3


        if np.isclose(ell_1,ell_2, atol = 0.00095) and np.isclose(ell_2, ell_3, atol  = 0.00095) and np.isclose(ell_1, ell_3, atol = 0.00095):
            theta = 0.0
            kappa = 0.0
            phi = 0

        else:
            kappa = 2*math.sqrt(ell_1**2 + ell_2**2 + ell_3**2 - ell_1 * ell_2 - ell_1*ell_3 - ell_2*ell_3)/(self.d*(ell_1 + ell_2 + ell_3))
            theta = kappa * S
            phi = math.atan2(math.sqrt(3) * (ell_2 + ell_3 - 2*ell_1), (3*(ell_2 - ell_3)))
            phi = phi % (2 * math.pi)

        for i in range(self.n):
            idx = i*3
            q_0[idx] = phi
            q_0[idx+1] = theta/self.n
            q_0[idx+2] = S/self.n

        return q_0

    def Compute_Jacobian_Ell_Q(self, q):

        # Function to evaluate the jacobian relating change in Ell w.r.t change in q

        q = np.asarray(q).reshape(-1)
        J = np.zeros((3,3*self.n))
        def Compute_d_Ell(q,d,i,j):
            if i ==0:
                dphi   = -d*q[j+1]*np.cos(q[j])
                dtheta = -d*np.sin(q[j-1])
            elif i ==1:
                dphi   = d*q[j+1]*np.sin(math.pi/6 - q[j])
                dtheta = d*np.cos(math.pi/6 - q[j-1])
            else:
                dphi   = d*q[j+1]*np.sin(math.pi/6 + q[j])
                dtheta = -d*np.cos(math.pi/6 + q[j-1])
            return dphi, dtheta

        for i in range(3):
            for j in range(0,3*self.n,3):
                J[i,j],_ = Compute_d_Ell(q,self.d,i,j)
            for j in range(1,3*self.n,3):
                _,J[i,j] = Compute_d_Ell(q,self.d,i,j)
            for j in range(2,3*self.n,3):
                J[i,j] = 1
        return J

    def Actuated_Jacobian(self, q):
        # Final actuated Jacobian equating change in tip_coordinates based on changes in Ell
        J_rq = self.ComputeJacobian(q)
        J_Eq = self.Compute_Jacobian_Ell_Q(q)

        J_Eq_inv = np.linalg.pinv(J_Eq)

        J_a = J_rq @ J_Eq_inv

        return J_a

    # Next three functions basically help to evaluate the plotting points based on the q, in other words
    # is the forward_kinematics of the SBA

    def Transform(self, theta, phi,l):
        T  = np.eye(4)

        H1 = np.eye(4)
        H1[:3, :3] = self.Rz(phi)

        H3 = np.eye(4)
        H3[:3, :3] = self.Rz(-phi)

        if abs(theta) < 1e-4:
            H2 = np.eye(4)
            H2[:3, :3] = self.Ry(theta)
            H2[:3,  3] = [0.0, 0.0, l]

            T = H1 @ H2 @ H3
        else:
            r = l/theta

            H2 = np.eye(4)
            H2[:3, :3] = self.Ry(theta)
            H2[:3,  3] = [r*(1 - np.cos(theta)),
                                            0.0,
                                r*np.sin(theta)]
            T = H1 @ H2 @ H3

        return T

    def Position_Vector(self, theta, samples, phi, l, T):
        t = np.linspace(0, theta, samples)
        X = np.zeros((3, samples))
        L = np.linspace(0,l,samples)

        for i in range(samples):
            Mat = self.Transform(t[i], phi, L[i])
            Transform_matrix = T @ Mat
            X[:, i] = Transform_matrix[:3, 3]

        return X

    def Compute_Soft_Curve (self, q):
        q = np.asarray(q).reshape(-1)
        samples = 100
        T = np.eye(4)
        softPts = np.zeros((3, samples*self.n))
        skelPts = np.zeros((3, self.n + 1))
        v = np.zeros((3,1))
        phi_s   = q[0::3]
        theta_s = q[1::3]
        L       = q[2::3]

        for i in range(self.n):
            Xi = self.Position_Vector(theta_s[i], samples, phi_s[i], L[i], T)
            start = i*samples
            end   = start + samples
            softPts[:, start:end] = Xi
            T = T @ self.Transform(theta_s[i], phi_s[i], L[i])
            skelPts[:, i+1] = softPts[:, end - 1]

        v = softPts[:, -1]

        return softPts, skelPts, v

    def Compute_actual_tip(self,q):
        dh_aug = self.Augmented_DH(q)
        new_dh_aug = dh_aug
        new_dh_aug[-1,1] = self.r

        T = np.eye(4)

        for k in range(len(new_dh_aug[:,0])):
            theta, d, a, alpha = new_dh_aug[k]
            T = T @ self.DH_transform(theta, d, a, alpha)

        actual_tip = T[:3,3]

        return actual_tip

    def ell_q(self,q):
        ell_1 = q.copy()
        ell_2 = q.copy()
        ell_3 = q.copy()
        for i in range(0,3*self.n,3):
            ell_1[2+i] = q[2+i] - self.d*q[i+1]*np.sin(q[i])
            ell_2[2+i] = q[2+i] + self.d*q[i+1]*np.cos((math.pi/6) - q[i])
            ell_3[2+i] = q[2+i] - self.d*q[i+1]*np.cos((math.pi/6) + q[i])

        return ell_1, ell_2, ell_3

    def ell(self,q):
        ell1 = np.sum(q[2::3] - self.d*q[1::3]*np.cos(math.pi/2 - q[::3]))
        ell2 = np.sum(q[2::3] - self.d*q[1::3]*np.cos(math.pi + (math.pi/6 - q[::3])))
        ell3 = np.sum(q[2::3] - self.d*q[1::3]*np.cos(2*math.pi - (math.pi/6 + q[::3])))

        return np.array([ell1, ell2, ell3])

    def force_change(self, q, del_uf):
        J_l = self.Actuated_Jacobian(q).T
        U, S, Vh = np.linalg.svd(J_l, full_matrices = False)
        V = Vh.T
        sigma_0 = 0.01*max(S)
        nu = 50
        h = (S**3 +nu*S**2 + 2*S + 2*sigma_0)/(S**2 + nu*S + 2)
        H_inv = np.diag(1.0/h)
        J_l_inv = V @ H_inv @ U.T

        delta_f = J_l_inv @ (0.6*del_uf)

        return delta_f

    def force_reading (self, force):
        rng = np.random.default_rng(seed=42)
        noise = rng.normal(loc=0.0, scale=self.sensor_var, size=np.shape(force))

        return force + noise

    def Jacobian_inverse(self, Jacob):
        #Findinf the inverse of non-square matrix using the Singular value filtering method
        U, S, Vh = np.linalg.svd(Jacob, full_matrices = False)
        V = Vh.T
        sigma_0 = 0.01*max(S)
        nu = 50
        h = (S**3 +nu*S**2 + 2*S + 2*sigma_0)/(S**2 + nu*S + 2)
        H_inv = np.diag(1.0/h)
        J_inv = V @ H_inv @ U.T

        return J_inv

def Plane_points(L, H, W, T_A):
    x_axis = [-L, L]
    points = np.zeros((8,4))
    c = 0
    for i_idx, i in enumerate(x_axis):
        for j_idx in range(len(x_axis)):
            points[c,:] = np.array([i,H*(-1)**((i_idx+1) + (j_idx + 1)),0,1]).reshape(1,-1)
            points[c+4,:] = np.array([i,H*(-1)**((i_idx+1) + (j_idx + 1)),W,1]).reshape(1,-1)
            c += 1
    vertices = T_A @ points.T
    vertices = vertices.T[:,:3]
    return vertices

def Path_line(t, Path, num_goal):
    Path_points = np.vstack([Path(t,j).ravel() for j in range(num_goal)])

    return Path_points

def isContact(tip_position, goal, plane_norm):
    if ((goal - tip_position).reshape(-1,1)).T@plane_norm.reshape(-1,1) <= 0:
        return True
    else:
        False


def InvKinematics(frame, plane_surf, error, vel, force_x, force_y, force_z, force_abs,
                  dist_plane, goal_mark, path_line, hSoft, Ell_1, Ell_2, Ell_3, hSkel,
                  top, rigid_portion, bottom, hForce,time_text, len1_mark,
                  len2_mark, len3_mark, kin, Path, Path_dot,Path_d_dot, Plane_origin):

    global i_goal
    global u_ells
    global time_start, time_prev
    global error_list, Time_list, Velocity_list, Force_list, Force_abs, dist
    global tip_no_load_prev, tip_model_prev
    global tip_model_vel
    global Ex_current, Ef_current
    global q_no_load_prev, q_model_prev
    global force_switch
    global num_goal, on_plane, stage
    global F_sensor, error_x, error_f, delta_uf_prev, plane_norm
    global F_desired, F_Measured
    global Projection_matrix, c, error_x_prev, delta_ux_prev, F_sensor_reading

    #To get the elapsed time and also the dt
    # time_now = time.time()
    # t_elapsed = time_now - time_start
    # dt = time_now - time_prev
    # time_prev = time_now

    # #time elapsed , dt is fixed #20MS
    dt = DT
    t_elapsed = frame * dt

    #Measure Force
    #F_sensor_reading = read_force_sensor()

    q_no_load_current = kin.q_no_load(u_ells)
    print(u_ells)
    q_no_load_dot     = kin.q_diff(q_no_load_current, q_no_load_prev)/dt

    q_model_current   = kin.q_model_dynamics(q_model_prev, q_no_load_current, q_no_load_dot, F_sensor, t_elapsed, dt)

    #q_model_corrected
    Plane_point = Path(t_elapsed, i_goal)
    ell_prev = kin.ell(q_model_prev).ravel()
    J_q = kin.ComputeJacobian(q_model_prev)
    J_ell = kin.Compute_Jacobian_Ell_Q(q_model_prev)

    F_ext = F_sensor

    print('F_ext: {}' .format(F_ext))

    A = kin.K + kin.D/dt
    b = kin.K @ q_no_load_current + kin.D @ q_no_load_dot + J_q.T @ F_ext + kin.D @ q_model_prev/dt

    ell_prev = np.asarray(ell_prev).ravel()
    u_ells   = np.asarray(u_ells).ravel()
    q_prev   = np.asarray(q_model_prev).ravel()
    q_curr   = np.asarray(q_model_current).ravel()
    esp = -0.05
    A_plane = (plane_norm.T @ J_q).reshape(1, -1)
    b_plane = float(plane_norm.T @ (Plane_point - tip_model_prev)
                -esp + (plane_norm.T @ J_q) @ q_prev)

    try:
        q = cp.Variable(q_prev.size)

        C_plane = (A_plane @ q <= b_plane)

        Cons = [
            plane_norm.T @ (Plane_point - (tip_model_prev + J_q @ (q - q_prev))) >= -0.05,
            # C_plane,
            (ell_prev - u_ells) + J_ell @ (q - q_prev) >= -0.2,
            (ell_prev - u_ells) + J_ell @ (q - q_prev) <= 0
        ]

        # Obj = cp.Minimize(0.5 * cp.sum_squares(q - q_model_current))

        beta = 1.0                      # <1 => more bending; >1 => stiffer correction
        W = beta * kin.K                # K is diagonal SPD; keeps correction stiffness-aware
        rho = 1e-2                       # dynamics penalty; lower => more bending, higher => closer to dynamics

        Obj = cp.Minimize(
            0.5 * cp.quad_form(q - q_model_current, W) +
            0.5 *cp.sum_squares(A@q - b)
        )

        q.value = q_curr
        problem = cp.Problem(Obj, Cons)
        problem.solve(solver=cp.OSQP, warm_start=True, eps_abs=1e-5, eps_rel=1e-5, max_iter=20000)




        if problem.status not in ("optimal", "optimal_inaccurate") or (q.value is None):
            print("QP status:", problem.status)
            q_model_corrected = q_curr.copy()

        else:
            q_model_corrected = q.value



    except Exception as e:
        print("Optimisation failed:", e)
        q_model_corrected = q_curr.copy()

    tau_i = (kin.K @ q_no_load_current + kin.D @ q_no_load_dot).ravel()
    tau_c = (A @ q_model_corrected) - tau_i - (kin.D @ q_model_prev)/dt
    # J = kin.ComputeJacobian(q_model_corrected)
    # reg = 1e-8
    # F_ext = np.linalg.solve(J @ J.T + reg*np.eye(3), J @ tau_c)

    J_ac = kin.ComputeJacobian(q_model_corrected)
    U, S, Vh = np.linalg.svd(J_ac, full_matrices = False)
    V = Vh.T
    sigma_0 = 0.01*max(S)
    nu = 50
    h = (S**3 +nu*S**2 + 2*S + 2*sigma_0)/(S**2 + nu*S + 2)
    H_inv = np.diag(1.0/h)
    J_ac_inv = V @ H_inv @ U.T
    F_ext = J_ac_inv.T @ tau_c

    gap_true = float(plane_norm.T @ (Plane_point - kin.Compute_actual_tip(q_model_corrected)))
    if gap_true >= 0.0:
        F_ext = np.zeros(3)
        # F_sensor = (-lam * plane_norm).ravel()
    F_sensor = F_ext

    print('--------------------------------------------------------------------------------------------------')
    print(c)
    print(i_goal)
    print('U_ells from correction: {}' .format(kin.ell(q_model_corrected)))
    print('U_ell given in correction: {}' .format(kin.ell(q_no_load_current)))
    print('U_ell without correction: {}' .format(kin.ell(q_model_current)))
    print('U_ell actual: {}' .format(u_ells))
    print('Testing: {}' .format(kin.ell(kin.q_no_load(u_ells))))

    print('q_no_load_current: {}' .format(q_no_load_current))
    print('q_model_current: {}' .format(q_model_current))
    print('q_model_corrected: {}' .format(q_model_corrected))
    print('Force from algo: {}' .format(F_ext))


    #try:
    #     def objective(q):
    #         dq = kin.q_diff(q, q_model_prev)/dt
    #         J = kin.ComputeJacobian(q)
    #         r = kin.K @ q + kin.D @ dq - b - J.T @ (F_ext)
    #         return r.dot(r)

    #     def dynamic_non_penetration(q):
    #         tip_pos = kin.Compute_actual_tip(q)
    #         return plane_norm.T @ (Plane_point - tip_pos)

    #     constraints = []
    #     constraints.append({
    #         'type': 'ineq',
    #         'fun': dynamic_non_penetration
    #     })
    #     q_corrected_guess = (q_model_current.ravel()).copy()

    #     res_q = minimize(
    #             objective,
    #             q_corrected_guess,
    #             method= 'SLSQP',
    #             constraints = constraints,
    #             options={'ftol':0.05, 'maxiter':20}
    #     )

    #     if not res_q.success:
    #         print("Solver Warning:", res_q.message)
    #         q_model_corrected= q_corrected_guess.copy()
    #     else:
    #         q_model_corrected = (res_q.x).ravel()

    # except Exception as e:
    #     print("Optimisation failed:", e)
        # q_model_corrected = q_corrected_guess.copy()

    # q_model_corrected = q_model_current.copy()

    tip_no_load_current = kin.Compute_actual_tip(q_no_load_current)
    tip_model_current   = kin.Compute_actual_tip(q_model_corrected)

    tip_model_vel       = (tip_model_current - tip_model_prev)/dt

    normal_distance     = plane_norm.T @ (Path(t_elapsed, i_goal) - tip_model_current)

#-----------------------------------------------------Plotting Portion ------------------------------------------------------------------
    plane_origin = Plane_origin(t_elapsed).reshape(-1,1)
    T_A_curr = np.block([[kin.Ry(kin.theta_plane), plane_origin],[np.zeros((1,3)), np.array([[1]])]])

    vertices = Plane_points(20,20,5,T_A_curr)

    faces = [
        [vertices[j] for j in [0,1,2,3]],  # bottom
        [vertices[j] for j in [4,5,6,7]],  # top
        [vertices[j] for j in [0,1,5,4]],  # front
        [vertices[j] for j in [1,2,6,5]],  # right
        [vertices[j] for j in [2,3,7,6]],  # back
        [vertices[j] for j in [3,0,4,7]]   # left
    ]

    plane_surf.set_verts(faces)

    ax = plane_surf.axes
    if ax.M is None:
        ax.get_proj()
    plane_surf.do_3d_projection()


    Path_way = lambda t: np.vstack([Path(t, j) for j in range(num_goal)])
    Path_list = Path_way(t_elapsed)

    if Path_list[-1,:].all != Path_list[0,:].all:
        Path_list = np.vstack((Path_list, Path_list[0,:].reshape(1,-1)))
    #Plotting of SBA and environment for current frame
    softPts, skelPts, _ = kin.Compute_Soft_Curve(q_model_corrected)

    ell_1, ell_2, ell_3 = kin.ell_q(q_model_corrected)
    ell_1_pts,_,_ = kin.Compute_Soft_Curve(ell_1)
    ell_2_pts,_,_ = kin.Compute_Soft_Curve(ell_2)
    ell_3_pts,_,_ = kin.Compute_Soft_Curve(ell_3)

    ell_1_origin = np.array([0,kin.d,0]).reshape(-1,1)
    ell_2_origin = np.array([-kin.d*np.cos(math.pi/6), -kin.d*np.sin(math.pi/6), 0]).reshape(-1,1)
    ell_3_origin = np.array([kin.d*np.cos(math.pi/6), -kin.d*np.sin(math.pi/6),0]).reshape(-1,1)


    ell_1_pts = ell_1_pts + ell_1_origin
    ell_2_pts = ell_2_pts + ell_2_origin
    ell_3_pts = ell_3_pts + ell_3_origin

    #Force vector from SBA to environment
    f_vec = -F_sensor

    #function handle for pathline/ pathcoords
    path_line.set_data(Path_list[:,0], Path_list[:,1])
    path_line.set_3d_properties(Path_list[:,2])

    #function handle for goal_point
    goal_mark.set_data([Path(t_elapsed, i_goal)[0]], [Path(t_elapsed, i_goal)[1]])
    goal_mark.set_3d_properties([Path(t_elapsed, i_goal)[2]])

    #function handle for soft backbone
    hSoft.set_data(softPts[0,:], softPts[1,:])
    hSoft.set_3d_properties(softPts[2,:])

    #function handle for chambers
    Ell_1.set_data(ell_1_pts[0,:], ell_1_pts[1,:])
    Ell_1.set_3d_properties(ell_1_pts[2,:])
    Ell_2.set_data(ell_2_pts[0,:], ell_2_pts[1,:])
    Ell_2.set_3d_properties(ell_2_pts[2,:])
    Ell_3.set_data(ell_3_pts[0,:], ell_3_pts[1,:])
    Ell_3.set_3d_properties(ell_3_pts[2,:])

    #function handle for RRPRR links
    hSkel.set_data(skelPts[0,:], skelPts[1,:])
    hSkel.set_3d_properties(skelPts[2,:])

    #function handle for top surface of SBA
    top.set_data([ell_1_pts[0,-1], ell_2_pts[0,-1], ell_3_pts[0,-1], ell_1_pts[0,-1]],
                 [ell_1_pts[1,-1], ell_2_pts[1,-1], ell_3_pts[1,-1], ell_1_pts[1,-1]])
    top.set_3d_properties([ell_1_pts[2,-1], ell_2_pts[2,-1], ell_3_pts[2,-1], ell_1_pts[2,-1]])

    #function handle for bottom surface of SBA
    bottom.set_data([ell_1_pts[0,0], ell_2_pts[0,0], ell_3_pts[0,0], ell_1_pts[0,0]],
                    [ell_1_pts[1,0], ell_2_pts[1,0], ell_3_pts[1,0], ell_1_pts[1,0]])
    bottom.set_3d_properties([ell_1_pts[2,0], ell_2_pts[2,0], ell_3_pts[2,0], ell_2_pts[2,0]])

    #function handle for rigid extended portion
    rigid_portion.set_data([skelPts[0,-1], tip_model_current[0]],
                           [skelPts[1,-1], tip_model_current[1]])
    rigid_portion.set_3d_properties([skelPts[2,-1], tip_model_current[2]])

    #function handle for force vector plot
    hForce.set_data([tip_model_current[0], tip_model_current[0] + 5*f_vec[0]],
                    [tip_model_current[1], tip_model_current[1] + 5*f_vec[1]])
    hForce.set_3d_properties([tip_model_current[2], tip_model_current[2] + 5*f_vec[2]])

    error.set_data(Time_list, error_list)
    vel.set_data(Time_list, Velocity_list)

    # function hanfle for chamber's length
    len1_mark.set_data([1,1], [0, u_ells[0]])
    len2_mark.set_data([2,2], [0, u_ells[1]])
    len3_mark.set_data([3,3], [0, u_ells[2]])

    force_x.set_data(Time_list, Force_list[:,0])
    force_y.set_data(Time_list, Force_list[:,1])
    force_z.set_data(Time_list, Force_list[:,2])

    force_abs.set_data(Time_list, Force_abs)

    dist_plane.set_data(Time_list, dist)

    #function handle for time
    time_text.set_text(f'Time: {t_elapsed:.2f}s')

#---------------------Control and modeling for next timeframe based on current spatial and force error---------------------
    # Parameters for motion control
    T = 0.1
    gamma = 0.1/T

    error_x = (Path(t_elapsed, i_goal) - tip_model_current)

    P_n = plane_norm @ plane_norm.T
    P_t = np.eye(3) - P_n

    distance  = float(plane_norm.T @ (Path(t_elapsed, i_goal) - tip_model_current))

    if (not on_plane) and distance < 0.4:
        on_plane = True
        Projection_matrix = P_n

    elif on_plane and distance > 0.6:
        on_plane = False
        Projection_matrix = np.zeros((3,3))


    #Measuring the force error between the sensor and the desired value
    error_f = Projection_matrix@(F_desired - F_sensor_reading)

    if np.linalg.norm(Path(t_elapsed,i_goal) - tip_model_current) < 0.4:
        i_goal = (i_goal + 1)%num_goal
        # i_goal=0
        # Ex_current = np.zeros(3)
        # error_x_prev = error_x

    # error_x = (np.eye(3) - Projection_matrix) @ error_x
    Ex_current += dt*error_x
    Ef_current += dt*error_f
    error_x_dot = (Path_dot(t_elapsed) - tip_model_vel)

    #Sum of previous error_f * dt
    # error_x_dot = (error_x - error_x_prev)/dt

    # print('Error_dot: {}' .format(error_x_dot))

    print('on_plane?: {}' .format(on_plane))
    print('X_error: {}' .format(error_x))
    print('F_error: {}' .format(error_f))


    #Getting the current Jacobian for the actuated coordinates ell
    J_ac = kin.Actuated_Jacobian(q_model_corrected)
    U, S, Vh = np.linalg.svd(J_ac, full_matrices = False)
    V = Vh.T
    sigma_0 = 0.01*max(S)
    nu = 50
    h = (S**3 +nu*S**2 + 2*S + 2*sigma_0)/(S**2 + nu*S + 2)
    H_inv = np.diag(1.0/h)
    J_ac_inv = V @ H_inv @ U.T

    # delta_f = Projection_matrix @ (0.8 * error_f + 0.1 * Ef_current)

    # delta_u = J_ac_inv @ (dt*Path_dot(t_elapsed) + dt*(14*error_x + 0.5*Ex_current + 0.8*error_x_dot)) + 0*J_ac.T @ (delta_f)

    J_ac = kin.Actuated_Jacobian(q_model_corrected)
    try:

        delta_u = cp.Variable(3)
        lam1 = 1
        delta_f = Projection_matrix @ (0.8 * error_f + 0.1 * Ef_current)

        constraints = [
            J_ac @ delta_u == (dt*Path_dot(t_elapsed) + dt*(13*error_x  + 0.09*error_x_dot))
        ]

        objective  = cp.Minimize(lam1*cp.norm(delta_u,2))

        prob = cp.Problem(objective, constraints)

        prob.solve()

        delta_u_val  = delta_u.value

    except (np.linalg.LinAlgError, cp.error.SolverError) as e:
        print("control optimisation failed:", e)
        delta_u_val = np.zeros(3)




    u_ells = (u_ells + kin.gain*delta_u_val).ravel()
    u_ells = np.round(np.clip(u_ells,2,50),3)

    #updating the ells based on the output from the control algorithm
    # del_f_measured = kin.force_change(q_current, delta_uf)
    # F_Measured = F_Measured  + del_f_measured
    # F_sensor_reading = (plane_norm@plane_norm.T) @F_Measured

    F_sensor_reading = np.zeros(3)

    tip_model_prev = tip_model_current.copy()
    tip_no_load_prev = tip_no_load_current.copy()
    q_model_prev = q_model_corrected.copy()
    q_no_load_prev = q_no_load_current.copy()

    error_x_prev = error_x


    Force_list = np.vstack((Force_list, -F_sensor))
    Force_abs.append([np.linalg.norm(-F_sensor)])
    error_list.append(np.linalg.norm(Path(t_elapsed, i_goal) - tip_model_current))
    Time_list.append(t_elapsed)
    Velocity_list.append(np.linalg.norm(tip_model_vel))
    dist.append([normal_distance])

    c+=1

    return plane_surf, error,vel,goal_mark, force_x, force_y, force_z, force_abs, dist_plane, path_line, hSoft, Ell_1, Ell_2, Ell_3, hSkel,top,rigid_portion, bottom, hForce, time_text, len1_mark, len2_mark, len3_mark

def sample_edge_points(V4x3, total):
    # xyz vertices
    A, B, C = V4x3[:3, :].T  # (3,) each

    # distribute counts per edge (AB, BC, CA)
    n1 = total // 3 + (total % 3 > 0)   # first edge
    n2 = total // 3 + (total % 3 > 1)   # second edge
    n3 = total - n1 - n2                # third edge

    def edge_linspace(P, Q, n):
        t = np.linspace(0, 1, n, endpoint=True)
        return P[None, :] + t[:, None] * (Q - P)[None, :]

    pts_ab = edge_linspace(A, B, n1)
    pts_bc = edge_linspace(B, C, n2)
    pts_ca = edge_linspace(C, A, n3)

    P = np.vstack([pts_ab, pts_bc, pts_ca])  # (50,3)
    # add homogeneous row
    P_h = np.vstack([P.T, np.ones(P.shape[0])])  # 4x50
    return P_h

def main():
    global u_ells, time_prev, goal_threshold, time_start
    global error_list, Time_list, Velocity_list
    global tip_no_load_prev, tip_model_prev
    global tip_model_vel
    global q_no_load_prev, q_model_prev
    global Force_list, Force_abs, dist
    global force_switch, error_f, error_x, plane_norm, F_sensor
    global num_goal, i_goal
    global F_desired
    global error_x_prev

    error_list =    []
    Time_list  =    []
    Velocity_list = []
    Force_abs     = []
    dist = []
    Force_list =    np.empty((0,3))


    kin = Kinematics()


    #Path points on the plane #w.r.t plane's origib
    Points_on_plane = np.array([[5, -5, -5],
                                [0,    -6.50,  6.50],
                                [0,        0,     0],
                                [1,        1,     1]])

    #defining the plane_normal
    plane_norm = kin.Ry(kin.theta_plane)[:3,2].reshape(-1,1)

    #Transformation matrix between plane and origin frame
    T_A_0 = np.block([[kin.Ry(kin.theta_plane), kin.P_plane_origin.reshape(-1,1)],[np.zeros((1,3)), np.array([[1]])]])

    #Path_coordinates in global reference
    Points_on_plane = sample_edge_points(Points_on_plane, 100)
    Path_Coords = T_A_0 @ Points_on_plane
    num_goal    = Path_Coords.shape[1]

    Plane_origin = lambda t: kin.P_plane_origin.ravel() + 5*plane_norm.ravel()*np.sin(math.pi*t/2)

    Path       = lambda t, i_goal: Path_Coords[:3, i_goal] + 5*plane_norm.ravel()*np.sin(math.pi*t/2)
    Path_dot   = lambda t: 5*plane_norm.ravel()*np.cos(math.pi*t/2)*math.pi/2
    Path_d_dot = lambda t: -5*plane_norm.ravel()*np.sin(math.pi*t/2)*(math.pi/2)**2


    #0.5N along the normal of the surface of interest
    f_given = np.array([0,0,1.0])
    F_given = kin.Ry(kin.theta_plane)@f_given

    F_desired = F_given

    #Length of each chambers at t = 0
    u_ells = np.array([10.3, 10.1, 10.4])

    #Initialization of shape of robot under no_load condition
    q_no_load_prev      = kin.q_no_load(u_ells)
    q_model_prev        = q_no_load_prev.copy()
    tip_no_load_prev    = kin.Compute_actual_tip(q_no_load_prev)
    tip_model_prev      = kin.Compute_actual_tip(q_model_prev)

    # At the start, there is no rate of change of config variables
    tip_model_vel   = 0.0
    tip_no_load_vel = 0.0

    #Normal distance from Surface
    Distance = plane_norm.T @ (Path(0, i_goal) - tip_model_prev)

    #Force sensor reading at t = 0
    f_measured    = np.zeros(3)

    F_sensor = f_measured


    #Initial position error and Force error
    error_x = Path(0, i_goal) - tip_model_prev
    error_x_prev = error_x
    error_f = F_desired - F_sensor


    Force_list = np.vstack((Force_list, -F_sensor))
    Force_abs.append([0])
    error_list.append(np.linalg.norm(kin.P_plane_origin - tip_model_prev))
    Time_list.append(0.0)
    Velocity_list.append(tip_model_vel)

    x = -1 if plane_norm.T @ (Path(0, i_goal) - tip_model_prev) < -0.1 else 1
    in_plane = x*np.linalg.norm(Path(0, i_goal) - tip_model_prev)

    normal_distance = plane_norm.T @(Path(0,i_goal) - tip_model_prev)
    dist.append([normal_distance])


    # Plotting and figure handles

    fig = plt.figure(figsize = (12,5))
    fig_2, (ax3_1, ax3_2) = plt.subplots(
        nrows=1, ncols=2,         # 1×2 grid
        figsize=(8, 4),           # overall figure size
    )



    ax1 = fig.add_subplot(1,2,1, projection ='3d')
    ax2 = fig.add_subplot(1,2,2)

    vertices0 = Plane_points(20,20,5, T_A_0)

    faces0 = [
        [vertices0[j] for j in [0,1,2,3]],  # bottom
        [vertices0[j] for j in [4,5,6,7]],  # top
        [vertices0[j] for j in [0,1,5,4]],  # front
        [vertices0[j] for j in [1,2,6,5]],  # right
        [vertices0[j] for j in [2,3,7,6]],  # back
        [vertices0[j] for j in [3,0,4,7]]   # left
    ]


    plane_surf = Poly3DCollection(faces0, color='black', alpha=0.15, animated = True)
    ax1.add_collection3d(plane_surf)

    fig.canvas.draw()

    error,     = ax3_1.plot([], [], '--', color = 'red', label = 'error_norm', linewidth = 2)
    vel,       = ax3_1.plot([], [], '--', color = 'green', label = 'tip_velocity', linewidth = 2)

    force_x,  = ax3_2.plot([], [], '--', color = 'blue', label = 'force_x', linewidth = 2)
    force_y,  = ax3_2.plot([], [], '--', color = 'red', label = 'force_y', linewidth = 2)
    force_z,  = ax3_2.plot([], [], '--', color = 'black', label = 'force_z', linewidth = 2)

    force_abs, = ax3_2.plot([], [], '--', color = 'green', label = 'force_norm', linewidth = 3)


    dist_plane, = ax3_2.plot([],[], '--',color = 'black', label = 'Normal distance from plane', linewidth = 1)


    goal_mark, = ax1.plot([], [], [], '.', color = 'tab:red', markersize = 4, label = 'goal_point')
    path_line, = ax1.plot([], [], [], '--', color = 'tab:blue', label = 'path_coords', markersize = 5)
    hSoft,     = ax1.plot([], [], [], '-', color = 'tab:orange', label = 'PCC', linewidth = 3, )
    Ell_1,     = ax1.plot([], [], [], '-', color = 'tab:blue', label = 'Ell_1', linewidth = 3, alpha = 0.5)
    Ell_2,     = ax1.plot([], [], [], '-', color = 'tab:green', label = 'Ell_2', linewidth = 3, alpha = 0.5)
    Ell_3,     = ax1.plot([], [], [], '-', color = 'tab:brown', label = 'Ell_3', linewidth = 3, alpha = 0.5)
    hSkel,     = ax1.plot([], [], [], '--o', color = 'black', label = 'Rigid_structure', linewidth = 1)

    top,       = ax1.plot([], [], [], '--', color = 'black', linewidth = 1, alpha = 0.5)

    rigid_portion, = ax1.plot([], [], [], '-', color = 'black', label = 'Offset', linewidth = 3)

    bottom,    = ax1.plot([], [], [], '--', color = 'black', linewidth = 1, alpha = 0.5)

    hForce, = ax1.plot([], [], [], '-', color='cyan', linewidth=3, label='Force Vector')

    time_text = ax1.text2D(0.02, 0.95, '', transform = ax1.transAxes,fontsize=12, bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    len1_mark, = ax2.plot([], [], '-o', color = 'tab:blue', label = 'chamber_1', linewidth = 3)
    len2_mark, = ax2.plot([], [], '-o', color = 'tab:green', label = 'chamber_2', linewidth = 3)
    len3_mark, = ax2.plot([], [], '-o', color = 'tab:brown', label = 'chamber_3', linewidth = 3)

    ax1.set_xlim([-35, 35])
    ax1.set_ylim([-35, 35])
    ax1.set_zlim([0,45])
    ax1.legend()
    ax1.set_xlabel('X-axis (mm)')
    ax1.set_ylabel('Y-axis (mm)')
    ax1.set_zlabel('Z-axis (mm)')
    ax1.grid(True)


    ax2.set_xlim([0, 4])
    ax2.set_ylim([0,50])
    ax2.legend()
    ax2.set_xlabel('Chamber  #')
    ax2.set_ylabel('chamber lengths (mm)')

    ax3_1.set_xlim([0, 100])
    ax3_1.set_ylim([0, 20])
    ax3_1.legend()
    ax3_1.set_xlabel('Time')
    ax3_1.set_ylabel('Value')

    ax3_2.set_xlim([-5, 100])
    ax3_2.set_ylim([-1, 2])
    ax3_2.legend()
    ax3_2.set_xlabel('Time')
    ax3_2.set_ylabel('Value (N)')

    time_start = time.time()

    time_prev= time_start - 0.05


    anim = animation.FuncAnimation(
        fig, InvKinematics, fargs = (plane_surf, error, vel, force_x, force_y, force_z, force_abs, dist_plane, goal_mark, path_line,hSoft,Ell_1, Ell_2, Ell_3, hSkel,top, rigid_portion, bottom, hForce,
                                     time_text, len1_mark, len2_mark, len3_mark, kin, Path, Path_dot, Path_d_dot, Plane_origin),
                                     cache_frame_data= False, interval = 10, blit = True
    )

    # Function to stop/play the simulation using spacebar
    paused = {'flag': False}
    def on_press(event):
        if event.key == ' ' or event.key == 'space':
            if paused['flag']:
                anim.event_source.start()
                paused['flag'] = False
            else:
                anim.event_source.stop()
                paused['flag'] = True


    fig.canvas.mpl_connect('key_press_event', on_press)

    plt.tight_layout()
    plt.show()



if __name__ == '__main__':
    main()

