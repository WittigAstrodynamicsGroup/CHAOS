# Code for mathematical functions related to quaternions and vector operations
"""
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Author: Kash Saddul
Institution: Astrodynamics Research Group, 
                University of Southampton
Development period: 2020-2024
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


quaternions.py



This file contains generic functions for quaternion operations, 
such as frame rotations and conversions.

Currently implemented functions:
    - inertial_to_body: Transforms a vector from inertial to body-fixed frame.
    - body_to_inertial: Transforms a vector from body-fixed to inertial frame.
    - project_vector: Projects a vector onto an axis (gets the component along an axis).
    - get_euler_angles: Extracts Euler angles from a given quaternion.
    - quaternion_multiplication: Multiplies two quaternions (Hamilton product).
    - hamilton_product: Alternative implementation of quaternion multiplication.
    - offset_angle: Calculates the angular offset between two frame rotations.
    - DCM_to_quat: Converts a direction cosine matrix (DCM) to a quaternion.
    - dcm_to_quaternion: Another implementation of DCM to quaternion conversion.
    - quaternion_to_DCM: Converts a quaternion to a DCM.
"""
import numpy as np


from scipy.spatial.transform import Rotation as R


def inertial_to_body(p, quat):
    """
    Transforms a vector from the inertial frame to the body-fixed frame using a quaternion.

    Args:
        p (numpy.ndarray): Input vector in the inertial frame (3D)
        quat (numpy.ndarray): Quaternion representing the rotation (4D)

    Returns:
        numpy.ndarray: Rotated vector in the body-fixed frame (3D)
    """

    #normalise the quaternion to unit length
    q = quat / np.linalg.norm(quat)

    #build the rotation matrix
    Amatrix = np.empty((3,3))
    Amatrix[0][0] = q[0]**2 - q[1]**2 - q[2]**2 + q[3]**2
    Amatrix[0][1] = 2 * (q[0]*q[1] + q[2] * q[3])
    Amatrix[0][2] = 2 * (q[0]*q[2] - q[1]*q[3])

    Amatrix[1][0] = 2 * (q[0]*q[1] - q[2] * q[3])
    Amatrix[1][1] = -q[0]**2 + q[1]**2 - q[2]**2 + q[3]**2
    Amatrix[1][2] = 2 * (q[1]*q[2] + q[0]*q[3])

    Amatrix[2][0] = 2 * (q[0]*q[2] + q[1]*q[3])
    Amatrix[2][1] = 2 * (q[1]*q[2] - q[0]*q[3])
    Amatrix[2][2] = -q[0]**2 - q[1]**2 + q[2]**2 + q[3]**2

    t_p = np.transpose(p)

    column_vectors = np.matmul(Amatrix, t_p)

    rotated_vector = np.transpose(column_vectors)
    
    return rotated_vector



# def project_vector(A, r):
#     #project a vector A into r, i.e the component of A parallel to r

#     #unitise the vector r,
#     r = r / np.linalg.norm(r)
#     A_parallel = np.dot(A, r)
    
#     return A_parallel

def body_to_inertial(p, q):
    """
    Transforms a vector from the body-fixed frame to the inertial frame using a quaternion.

    This function performs the inverse of the `inertial_to_body` function. 
    It achieves this by applying the rotation represented by the inverse of the 
    provided quaternion (q_inv) to the vector (p).

    Args:
        p (numpy.ndarray): Input vector in the body-fixed frame (3D)
        q (numpy.ndarray): Quaternion representing the rotation (4D)

    Returns:
        numpy.ndarray: Rotated vector in the inertial frame (3D)
    """

    #inverse the Euler axis
    q_inv = np.array([-q[0], -q[1], -q[2], q[3]])
    
    #Apply same rotation, with inverse quaternion
    p_in = inertial_to_body(p, q_inv)              #Quaternion is normaised inside this function

    return p_in


def get_euler_angles( q):
    """
    Extracts Euler angles (precession, nutation, spin) from a given quaternion.

    This function calculates the Euler angles (alpha, beta, gamma) representing 
    the same rotation as the provided quaternion (q) using the direction cosine matrix (A).

    Args:
        q (numpy.ndarray): Quaternion representing a rotation (4D)

    Returns:
        numpy.ndarray: Euler angles (alpha, beta, gamma) in degrees (3D)
    """




    #direction cosine matrix A:
    Amatrix = np.empty((3,3))
    Amatrix[0][0] = q[0]**2 - q[1]**2 - q[2]**2 + q[3]**2
    Amatrix[0][1] = 2 * (q[0]*q[1] + q[2] * q[3])
    Amatrix[0][2] = 2 * (q[0]*q[2] - q[1]*q[3])

    Amatrix[1][0] = 2 * (q[0]*q[1] - q[2] * q[3])
    Amatrix[1][1] = -q[0]**2 + q[1]**2 - q[2]**2 + q[3]**2
    Amatrix[1][2] = 2 * (q[1]*q[2] + q[0]*q[3])

    Amatrix[2][0] = 2 * (q[0]*q[2] + q[1]*q[3])
    Amatrix[2][1] = 2 * (q[1]*q[2] - q[0]*q[3])
    Amatrix[2][2] = -q[0]**2 - q[1]**2 + q[2]**2 + q[3]**2


    #euler angles:
    alpha = np.arctan(Amatrix[2][0] / (-Amatrix[2][1]) ) % 2 * np.pi #precession

    beta = np.arccos(Amatrix[2][2]) % np.pi # nutation

    gamma = np.arctan(Amatrix[0][2] / Amatrix[1][2]) % 2 * np.pi # spin

    return np.array([alpha*180/np.pi, beta*180/np.pi, gamma*180/np.pi])         #in deg


def quaternion_multiplication(q1, q2):
    """
    Multiplies two quaternions using the Hamiltonian product.

    This function performs the quaternion multiplication (q1 * q2) according to the 
    Hamilton product, which involves separate calculations for the vector and scalar parts.

    Args:
        q1 (numpy.ndarray): First quaternion (4D)
        q2 (numpy.ndarray): Second quaternion (4D)

    Returns:
        numpy.ndarray: The product quaternion (4D)
    """
    #crassidis defines this product as - np.cross!! 
    #But real Hamiltonian product is +np.cross
    vec_part = q2[3] * q1[:3] + q1[3] * q2[:3] + np.cross(q1[:3], q2[:3])
    scalar_part = q1[3] * q2[3] - np.dot(q1[:3], q2[:3])

    quat = np.append(vec_part, scalar_part)

    return quat

def hamilton_product(q_a, q_b):
    """
    Calculates the Hamiltonian product of two quaternions.

    This function performs the element-wise multiplication of two quaternions (q_a and q_b) 
    following the rules of the Hamiltonian product. It separates the calculations for the 
    scalar and vector parts of the resulting product quaternion.

    Args:
        q_a (numpy.ndarray): First quaternion (4D) - denoted as a for clarity
        q_b (numpy.ndarray): Second quaternion (4D) - denoted as b for clarity

    Returns:
        numpy.ndarray: The product quaternion (4D)
    """
    a1, a2, a3, a4 = q_a
    b1, b2, b3, b4 = q_b

    scalar = a4*b4 - a1*b1 - a2*b2 - a3*b3

    vector1 = a1*b4 + a2*b3 - a3*b2 + a4*b1
    vector2 = -a1*b3 + a2*b4 + a3*b1 + a4*b2
    vector3 = a1*b2 - a2*b1 + a3*b4 + a4*b3

    prod = np.array([vector1, vector2, vector3, scalar])

    return prod



def offset_angle(q_ref1, q_in1):
    '''
    Calculates the in-plane angular offset between two rotations represented by quaternions.

    This function takes two quaternions (q_ref1 and q_in1) representing frame rotations 
    and computes the maximum angular offset between them. 

    '''
    #find error quaternion
    #inverse quaternion
    qx, qy, qz, qr = q_in1
    q_in1_inv = np.array([-qx, -qy, -qz, qr])

    q_err1 = quaternion_multiplication(q_ref1, q_in1_inv)
    q_err = q_err1/ np.linalg.norm(q_err1)
    q_err_norm = np.linalg.norm(q_err)


    #Principal angle of the error quaternion
    omega_err = 2 *np.arccos(q_err[3])
    err_axis = q_err[:3] / (np.sin(omega_err / 2))

    return omega_err


def DCM_to_quat(A):
    """
    Converts a direction cosine matrix (DCM) to a quaternion.

    This function takes a 3x3 direction cosine matrix (A) as input and computes the 
    corresponding quaternion representation. The quaternion is calculated using a 
    formula that exploits the trace of the matrix.

    Args:
        A (numpy.ndarray): Direction cosine matrix (3x3)

    Returns:
        numpy.ndarray: Quaternion representation (4D)
    """

    qx = np.sqrt(0.25 * ( 1 + A[0][0] - A[1][1] - A[2][2]))
    qy = np.sqrt(0.25 * ( 1 - A[0][0] + A[1][1] - A[2][2]))
    qz = np.sqrt(0.25 * ( 1 - A[0][0] - A[1][1] + A[2][2]))
    qr = np.sqrt(0.25 * ( 1 + A[0][0] + A[1][1] + A[2][2]))

    return np.array([qx, qy, qz, qr])



def dcm_to_quaternion(dcm):
    """
    Converts a direction cosine matrix (DCM) to a quaternion using eigenvalue decomposition.

    This function takes a 3x3 direction cosine matrix (dcm) as input and computes the 
    corresponding quaternion representation. It utilizes eigenvalue decomposition to find 
    the eigenvector corresponding to the largest eigenvalue, which is then used to construct 
    the quaternion. Finally, the quaternion is normalized to unit length. 

    Args:
        dcm (numpy.ndarray): Direction cosine matrix (3x3)

    Returns:
        numpy.ndarray: Quaternion representation (4D)
    """


    k = (1/3) * np.array([[dcm[0, 0] - dcm[1, 1] - dcm[2, 2], dcm[1, 0] + dcm[0, 1], dcm[2, 0] + dcm[0, 2], dcm[1, 2] - dcm[2, 1]],
                  [dcm[1, 0] + dcm[0, 1], dcm[1, 1] - dcm[0, 0] - dcm[2, 2], dcm[2, 1] + dcm[1, 2], dcm[2, 0] - dcm[0, 2]],
                  [dcm[2, 0] + dcm[0, 2], dcm[2, 1] + dcm[1, 2], dcm[2, 2] - dcm[0, 0] - dcm[1, 1], dcm[0, 1] - dcm[1, 0]],
                  [dcm[1, 2] - dcm[2, 1], dcm[2, 0] - dcm[0, 2], dcm[0, 1] - dcm[1, 0], dcm[0, 0] + dcm[1, 1] + dcm[2, 2]]])
    
    eigenvalues, eigenvectors = np.linalg.eig(k)
    
    # Find the eigenvector corresponding to the largest eigenvalue
    max_index = np.argmax(eigenvalues)
    quaternion = eigenvectors[:, max_index]
    
    # Normalize the quaternion to unit length
    quaternion /= np.linalg.norm(quaternion)
    
    return quaternion



def quaternion_to_DCM(quat):
    """
    Convert a quaternion to a 3x3 rotation matrix.
    
    Args:
        quat (list or numpy array): Quaternion [qx, qy, qz, qr]
    
    Returns:
        numpy array: 3x3 rotation matrix
    """
    qx, qy, qz, qr = quat

    # Compute quaternion normalization factor
    norm = np.sqrt(qx**2 + qy**2 + qz**2 + qr**2)
    if norm == 0:
        return np.eye(3)

    # Compute squared quaternion elements
    qx2 = qx**2
    qy2 = qy**2
    qz2 = qz**2
    qr2 = qr**2

    # Compute cross-products
    qx_qy = qx * qy
    qx_qz = qx * qz
    qx_qr = qx * qr
    qy_qz = qy * qz
    qy_qr = qy * qr
    qz_qr = qz * qr

    # Compute rotation matrix elements
    m11 = qx2 - qy2 - qz2 + qr2
    m12 = 2 * (qx_qy - qz_qr)
    m13 = 2 * (qx_qz + qy_qr)
    m21 = 2 * (qx_qy + qz_qr)
    m22 = -qx2 + qy2 - qz2 + qr2
    m23 = 2 * (qy_qz - qx_qr)
    m31 = 2 * (qx_qz - qy_qr)
    m32 = 2 * (qy_qz + qx_qr)
    m33 = -qx2 - qy2 + qz2 + qr2

    # Construct the rotation matrix
    rotation_matrix = np.array([[m11, m12, m13],
                                [m21, m22, m23],
                                [m31, m32, m33]])

    # Normalize the rotation matrix
    rotation_matrix /= norm

    return rotation_matrix

