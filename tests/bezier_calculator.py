import numpy as np

def cubic_bezier_curve(t, P0, P1, P2, P3):
	return (1-t)**3 * P0 + \
		   3 * (1-t)**2 * t * P1 + \
		   3*(1-t) * t**2 * P2 + \
		   t**3 * P3

def d_cubic_bezier_curve(t, P0, P1, P2, P3):
	return 3 * (1-t)**2 * (P1-P0) + \
		   6 * (1-t) * t * (P2-P1) + \
		   3 * t**2 * (P3-P2)

def scaled_d_cubic_bezier_curve(t, T, P0, P1, P2, P3):
	return 1/T * d_cubic_bezier_curve(t, P0, P1, P2, P3)

def d_d_cubic_bezier_curve(t, P0, P1, P2, P3):
	return 6 * (1-t) * (P2 - 2*P1 + P0) + \
		   6 * t * (P3 - 2*P2 + P1)

def scaled_d_d_cubic_bezier_curve(t, T, P0, P1, P2, P3):
	return (1/T)**2 * d_d_cubic_bezier_curve(t, P0, P1, P2, P3)

def quintic_bezier_curve(t, P0, P1, P2, P3, P4, P5):
	return (1-t)**5 * P0 + \
		   5*(1-t)**4 * t * P1 + \
		   10 * t**2 * (1-t)**3 * P2 + \
		   10 * t**3 * (1-t)**2 * P3 + \
		   5 * t**4 * (1-t) * P4 + \
		   t**5 * P5

####################### evaluate cubic #######################

def evaluate_cubic():
	P0 = np.array([10, 5])
	P1 = np.array([8, 2])
	P2 = np.array([6, 3])
	P3 = np.array([1, 20])
	print("cubic")
	for i in np.linspace(0, 1, 11):
		print("t: {}, value: {}".format(i, cubic_bezier_curve(i, P0, P1, P2, P3)))

def evaluate_d_cubic():
	P0 = np.array([10, 5])
	P1 = np.array([8, 2])
	P2 = np.array([6, 3])
	P3 = np.array([1, 20])
	print("d_cubic")
	for i in np.linspace(0, 1, 11):
		print("t: {}, value: {}".format(i, d_cubic_bezier_curve(i, P0, P1, P2, P3)))

def evaluate_scaled_d_cubic():
	P0 = np.array([10, 5])
	P1 = np.array([8, 2])
	P2 = np.array([6, 3])
	P3 = np.array([1, 20])
	print("scaled_d_cubic")
	T = 10
	for i in np.linspace(0, 1, 11):
		print("t: {}, value: {}".format(i, scaled_d_cubic_bezier_curve(i, T, P0, P1, P2, P3)))

def evaluate_d_d_cubic():
	P0 = np.array([10, 5])
	P1 = np.array([8, 2])
	P2 = np.array([6, 3])
	P3 = np.array([1, 20])
	print("d_d_cubic")
	for i in np.linspace(0, 1, 11):
		print("t: {}, value: {}".format(i, d_d_cubic_bezier_curve(i, P0, P1, P2, P3)))

def evaluate_scaled_d_d_cubic():
	P0 = np.array([10, 5])
	P1 = np.array([8, 2])
	P2 = np.array([6, 3])
	P3 = np.array([1, 20])
	print("scaled_d_d_cubic")
	T = 10
	for i in np.linspace(0, 1, 11):
		print("t: {}, value: {}".format(i, scaled_d_d_cubic_bezier_curve(i, T, P0, P1, P2, P3)))

####################### evaluate quintic #######################

def evaluate_quintic():
	P0 = np.array([15, 43])
	P1 = np.array([1, -5])
	P2 = np.array([10, -10])
	P3 = np.array([50, 32])
	P4 = np.array([-79, 23])
	P5 = np.array([5, -8])
	print("quintic")
	for i in np.linspace(0, 1, 11):
		print("t: {}, value: {}".format(i, quintic_bezier_curve(i, P0, P1, P2, P3, P4, P5)))

# evaluate_cubic()
# evaluate_d_cubic()
# evaluate_scaled_d_cubic()
# evaluate_d_d_cubic()
evaluate_scaled_d_d_cubic()
# evaluate_quintic()