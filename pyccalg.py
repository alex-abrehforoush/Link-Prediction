import os, sys, getopt, time, random
from math import log
import pulp as plp
from scipy.optimize import linprog

separator = '---------------'
eps = 0.0000000001

def _running_time_ms(start):
	return int(round((time.time()-start)*1000))

def _load(dataset_path,random_edgeweight_generation,edge_addition_prob):
	with open(dataset_path) as f:
		tot_min = 0
		id2vertex = {}
		vertex2id = {}
		edges = []
		graph = {}
		vertex_id = 0
		for line in f.readlines()[1:]:
			tokens = line.split()
			u = int(tokens[0])
			v = int(tokens[1])
			if random_edgeweight_generation:
				wp = random.uniform(random_edgeweight_generation[0],random_edgeweight_generation[1])
				wn = random.uniform(random_edgeweight_generation[0],random_edgeweight_generation[1])
			else:
				wp = float(tokens[2])
				wn = float(tokens[3])
				#wp = 1.0
				#wn = 0.0
			if wp != wn:
				if u not in vertex2id:
					vertex2id[u] = vertex_id
					id2vertex[vertex_id] = u
					vertex_id += 1
				if v not in vertex2id:
					vertex2id[v] = vertex_id
					id2vertex[vertex_id] = v
					vertex_id += 1
				uid = vertex2id[u]
				vid = vertex2id[v]
				if uid < vid:
					edges.append((uid,vid))
				else:
					edges.append((vid,uid))
				if uid not in graph.keys():
					graph[uid] = {}
				if vid not in graph.keys():
					graph[vid] = {}
				min_pn = min(wp,wn)
				tot_min += min_pn
				graph[uid][vid] = (wp-min_pn,wn-min_pn)
				graph[vid][uid] = (wp-min_pn,wn-min_pn)

		if random_edgeweight_generation and edge_addition_prob > 0:
			#generate weights for non-existing edges (with probability 'edge_addition_prob')
			for uid in graph.keys():
				for vid in id2vertex.keys():
					if uid<vid and vid not in graph[u] and random.random()<=edge_addition_prob:
						wp = random.uniform(random_edgeweight_generation[0],random_edgeweight_generation[1])
						wn = random.uniform(random_edgeweight_generation[0],random_edgeweight_generation[1])
						min_pn = min(wp,wn)
						tot_min += min_pn
						graph[uid][vid] = (wp-min_pn,wn-min_pn)
						graph[vid][uid] = (wp-min_pn,wn-min_pn)
						edges.append((uid,vid))

		return (id2vertex,vertex2id,edges,graph,tot_min)

def _read_params():
	dataset_file = None
	random_edgeweight_generation = None
	solver = 'scipy'
	algorithm = 'charikar'
	edge_addition_prob = -1
	short_params = 'd:r:s:a:m:'
	long_params = ['dataset=','random=','solver=','addedges=','method=']
	try:
		arguments, values = getopt.getopt(sys.argv[1:], short_params, long_params)
	except getopt.error as err:
		print('ologncc.py -d <dataset_file> [-r <rnd_edge_weight_LB,rnd_edge_weight_UB>] [-a <edge_addition_probability>] [-s <solver>] [-m <algorithm>]')
		sys.exit(2)
	for arg, value in arguments:
		if arg in ('-d', '--dataset'):
			dataset_file = value
		elif arg in ('-r', '--random'):
			random_edgeweight_generation = [float(x) for x in value.split(',')]
		elif arg in ('-s', '--solver'):
			solver = value.lower()
		elif arg in ('-a', '--addedges'):
			edge_addition_prob = float(value)
		elif arg in ('-m', '--method'):
			algorithm = value.lower()
	return (dataset_file,random_edgeweight_generation,edge_addition_prob,solver,algorithm)

def _map_cluster(cluster,id2vertex):
	return {id2vertex[u] for u in cluster}

def _vertex_pair_id(i,j,n):
	if i == j:
		raise Exception('ERROR: i and j must be different')
	lb = min(i,j)
	ub = max(i,j)
	return int(lb*n-(lb*(lb+1)/2))+ub-lb-1

def _vertex_pair_ids(n):
	id2vertexpair = {}
	id = 0
	for i in range(n-1):
		for j in range(i+1,n):
			id2vertexpair[id] = (i,j)
			id += 1
	return id2vertexpair

def _linear_program_scipy(num_vertices,edges,graph):
	vertex_pairs = int(num_vertices*(num_vertices-1)/2)
	A = []
	"""
	#old way of generating triangle-inequality constraints (deprecated as more elaborated and confusing, though still correct)
	for i in range(num_vertices-2):
		for j in range(i+1,num_vertices-1):
			for k in range(j+1,num_vertices):
				ik = _vertex_pair_id(i,k,num_vertices)
				ij = _vertex_pair_id(i,j,num_vertices)
				jk = _vertex_pair_id(j,k,num_vertices)
				#for all i,j,k, 3 triangle-inequality constraints should be stated:
				#First triangle-inequality constraint: xik <= xij + xjk <=> xik - xij - xjk = 0
				a = [0]*vertex_pairs
				a[ik] = 1
				a[ij] = -1
				a[jk] = -1
				A.append(a)
				#Second triangle-inequality constraint: xij <= xik + xjk <=> xij - xik - xjk = 0
				a = [0]*vertex_pairs
				a[ij] = 1
				a[ik] = -1
				a[jk] = -1
				A.append(a)
				#Third triangle-inequality constraint: xjk <= xij + xik <=> xjk - xij - xik = 0
				a = [0]*vertex_pairs
				a[jk] = 1
				a[ij] = -1
				a[ik] = -1
				A.append(a)
	"""
	for i in range(num_vertices-1):
		for j in range(i+1,num_vertices):
			ij = _vertex_pair_id(i,j,num_vertices)
			for k in range(num_vertices):
				if k != i and k != j:
					ik = _vertex_pair_id(i,k,num_vertices)
					kj = _vertex_pair_id(k,j,num_vertices)
					#for all vertex pairs {i,j} and all vertices k \notin {i,j}, state the following triangle-inequality constraint:
					# xij <= xik + xkj <=> xij - xik - xkj = 0
					a = [0]*vertex_pairs
					a[ij] = 1
					a[ik] = -1
					a[kj] = -1
					A.append(a)
	b = [0]*len(A)
	c = [0]*vertex_pairs
	for (u,v) in edges:
		uv = _vertex_pair_id(u,v,num_vertices)
		(wp,wn) = graph[u][v]
		if wp != wn:
			if wp < wn: #(u,v) \in E^-
				c[uv] = -(wn-wp)
			else: #(u,v) \in E^+
				c[uv] = (wp-wn)
	return (A,b,c)

def _solve_lp_scipy(A,b,c):
	#notes on supported solvers (see https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.linprog.html)
	#Method ‘highs-ds’ is a wrapper of the C++ high performance dual revised simplex implementation (HSOL) [13], [14].
	#Method ‘highs-ipm’ is a wrapper of a C++ implementation of an interior-point method [13];
	#it features a crossover routine, so it is as accurate as a simplex solver.
	#Method ‘highs’ chooses between the two automatically.
	#For new code involving linprog, we recommend explicitly choosing one of these three method values
	#instead of ‘interior-point’ (default), ‘revised simplex’, and ‘simplex’ (legacy).
	#lp_solution = linprog(c, A_ub=A, b_ub=b, bounds=[(0,1)])
	#lp_solution = linprog(c, A_ub=A, b_ub=b, bounds=[(0,1)], method='simplex')
	lp_solution = linprog(c, A_ub=A, b_ub=b, bounds=[(0,1)], method='highs-ipm')
	#lp_solution = linprog(c, A_ub=A, b_ub=b, bounds=[(0,1)], method='highs')
	lp_var_assignment = lp_solution['x']
	obj_value = lp_solution['fun']
	for i in range(len(lp_var_assignment)):
		if lp_var_assignment[i] < 0:
			lp_var_assignment[i] = 0
		elif lp_var_assignment[i] > 1:
			lp_var_assignment[i] = 1
	return (lp_var_assignment,obj_value)

def _linear_program_pulp(num_vertices,edges,graph):
	#see https://medium.com/opex-analytics/optimization-modeling-in-python-pulp-gurobi-and-cplex-83a62129807a
	opt_model = plp.LpProblem(name='GeneralWeightedCC')
	vertex_pairs = int(num_vertices*(num_vertices-1)/2)

	x_vars  = {i: plp.LpVariable(cat=plp.LpContinuous, lowBound=0, upBound=1, name='x_{0}'.format(i)) for i in range(vertex_pairs)}

	c_count = 0
	constraints = {}
	"""
	#old way of generating triangle-inequality constraints (deprecated as more elaborated and confusing, though still correct)
	for i in range(num_vertices-2):
		for j in range(i+1,num_vertices-1):
			for k in range(j+1,num_vertices):
				ik = _vertex_pair_id(i,k,num_vertices)
				ij = _vertex_pair_id(i,j,num_vertices)
				jk = _vertex_pair_id(j,k,num_vertices)
				#for all i,j,k, 3 triangle-inequality constraints should be stated:
				expr1 = plp.LpAffineExpression([(x_vars[ik],1), (x_vars[ij],-1), (x_vars[jk],-1)])#First triangle-inequality constraint: xik <= xij + xjk <=> xik - xij - xjk = 0
				expr2 = plp.LpAffineExpression([(x_vars[ij],1), (x_vars[ik],-1), (x_vars[jk],-1)])#Second triangle-inequality constraint: xij <= xik + xjk <=> xij - xik - xjk = 0
				expr3 = plp.LpAffineExpression([(x_vars[jk],1), (x_vars[ij],-1), (x_vars[ik],-1)])#Third triangle-inequality constraint: xjk <= xij + xik <=> xjk - xij - xik = 0
				for expr in [expr1,expr2,expr3]:
					constraints[c_count] = opt_model.addConstraint(plp.LpConstraint(e=expr,sense=plp.LpConstraintLE,rhs=0,name='constraint_{0}'.format(c_count)))
					c_count += 1
	"""
	for i in range(num_vertices-1):
		for j in range(i+1,num_vertices):
			ij = _vertex_pair_id(i,j,num_vertices)
			for k in range(num_vertices):
				if k != i and k != j:
					ik = _vertex_pair_id(i,k,num_vertices)
					kj = _vertex_pair_id(k,j,num_vertices)
					#for all vertex pairs {i,j} and all vertices k \notin {i,j}, state the following triangle-inequality constraint:
					# xij <= xik + xkj <=> xij - xik - xkj = 0
					expr = plp.LpAffineExpression([(x_vars[ij],1), (x_vars[ik],-1), (x_vars[kj],-1)])
					constraint = plp.LpConstraint(e=expr,sense=plp.LpConstraintLE,rhs=0,name='constraint_{0}'.format(c_count))
					constraints[c_count] = constraint
					opt_model.addConstraint(constraint)
					c_count += 1

	obj_expr = []
	for (u,v) in edges:
		uv = _vertex_pair_id(u,v,num_vertices)
		(wp,wn) = graph[u][v]
		if wp != wn:
			if wp < wn: #(u,v) \in E^-
				obj_expr.append(-(wn-wp)*x_vars[uv])
			else: #(u,v) \in E^+
				obj_expr.append((wp-wn)*x_vars[uv])
	objective = plp.lpSum(obj_expr)

	opt_model.sense = plp.LpMinimize
	opt_model.setObjective(objective)

	"""
	#######################
	#######################
	#######################
	## DEBUG:
	(A,b,c) = _linear_program_scipy(num_vertices,edges,graph)
	checkA = _check_constraint_correspondence(A,constraints)
	if not checkA:
		raise Exception('No constraint correspondence')

	checkc = _check_objective_correspondence(c,objective)
	if not checkc:
		raise Exception('No objective correspondence')
	#######################
	#######################
	#######################
	"""

	return opt_model

def _solve_lp_pulp(model):
	#model.solve()
	#model.solve(solver=plp.PULP_CBC_CMD(fracGap=0.00001))
	model.solve(solver=plp.PULP_CBC_CMD(msg=False))
	#lp_var_assignment = [x.varValue for x in model.variables()]
	lp_var_assignment = [0]*len(model.variables())
	for var in model.variables():
		varname = var.name
		varindex = int(varname.split('_')[1])
		lp_var_assignment[varindex] = var.varValue
	obj_value = model.objective.value()
	for i in range(len(lp_var_assignment)):
		if lp_var_assignment[i] < 0:
			lp_var_assignment[i] = 0
		elif lp_var_assignment[i] > 1:
			lp_var_assignment[i] = 1
	return (lp_var_assignment,obj_value)

def _sorted_distances(u,valid_vertices,num_vertices,x):
	du = []
	for v in valid_vertices:
		if v != u:
			du.append((v,x[_vertex_pair_id(u,v,num_vertices)]))
	return sorted(du,key=lambda f: f[1])

def _cut(ball,valid_vertices,graph):
	cut = 0
	for u in ball:
		if u in graph:
			for v in graph[u]:
				if v in valid_vertices and v not in ball:
					cut += graph[u][v][0]
	return cut

def _incremental_cut(old_ball,new_ball_vertices,valid_vertices,graph):
	incr_cut = 0
	for u in new_ball_vertices:
		if u in graph:
			for v in graph[u]:
				if v in valid_vertices:
					if v in old_ball:
						incr_cut -= graph[u][v][0]
					elif v not in new_ball_vertices:
						incr_cut += graph[u][v][0]
	return incr_cut

# need to pass the ball center ('u' in the paper) in order to properly compute the fractional weighted distance of positive edges leaving the ball
def _vol(ball_center,ball,valid_vertices,graph,num_vertices,x,r):
	vol = 0
	for v in ball:
		if v in graph:
			xuv = 0 if v==ball_center else x[_vertex_pair_id(ball_center,v,num_vertices)]
			#print('u: %d, v: %d, xuv: %s' %(ball_center,v,xuv))
			for w in graph[v]:
				if w in valid_vertices:
					idvw = _vertex_pair_id(v,w,num_vertices)
					xvw = x[idvw]
					if w in ball:
						if v<w: #check not to consider xvw two times
							vol += xvw*graph[v][w][0]
					else:
						if r < xuv:
							raise Exception('ERROR: radius cannot be less than the distance between the ball center and any vertex in the ball---r: %s, xuv: %s' %(r,xuv))
						vol += (r-xuv)*graph[v][w][0]
	return vol

#'F' in the paper
def _vol_whole_graph(graph,num_vertices,x):
	vol = 0
	for u in graph.keys():
		for v in graph[u]:
			if u < v: #check not to consider xuv two times
				xuv = x[_vertex_pair_id(u,v,num_vertices)]
				cuv = graph[u][v][0]
				vol += xuv*cuv
	return vol

#rounding algorithm proposed in Demaine et al., "Correlation clustering in general weighted graphs", TCS 2006
def round_demaine(x,id2vertexpair,id2vertex,edges,graph,const):
	clusters = []
	n = len(id2vertex)
	remaining_vertices = set(id2vertex.keys())
	shuffled_vertices = list(id2vertex.keys())
	random.shuffle(shuffled_vertices)
	F = _vol_whole_graph(graph,n,x)

	for u in shuffled_vertices:
		if u in remaining_vertices:
			du = _sorted_distances(u,remaining_vertices,n,x)
			initial_vol = F/n #default initial volume of a ball; while it is defined in the paper, it is just a technicality; the herein implementation allows for avoiding it

			#initial ball is composed of all vertices at distance 0 from u
			ball = {u}
			i = 0
			while i<len(du) and du[i][1]==0:
				i += 1
			for (v,d) in du[0:i]:
				ball.add(v)
			cut = _cut(ball,remaining_vertices,graph)
			#vol = initial_vol + _vol(u,ball,remaining_vertices,graph,n,x,1/const)
			du = du[i:]

			while du:
				r = min(1/const,du[0][1]) #du[0][1]: minimum distance in du
				i = 0
				while i<len(du) and du[i][1]<=r:
					i += 1
				new_ball_vertices = {v for (v,d) in du[0:i]}
				if new_ball_vertices:
					#'cut', 'ball', and 'du' are updated only if new ball vertices have been found
					cut += _incremental_cut(ball,new_ball_vertices,remaining_vertices,graph) #cut can be computed incrementally
					ball.update(new_ball_vertices)
					du = du[i:]
				#'vol' is updated anyway, as 'r' may have changed
				#it is not convenient to compute 'vol' incrementally as r changes in every iteration, so the vertices in the current ball must be all visited again anyway
				vol = initial_vol + _vol(u,ball,remaining_vertices,graph,n,x,r)

				if not new_ball_vertices or cut <= const*log(n+1)*vol: #'log' returns natural logarithm
					if not new_ball_vertices and cut > const*log(n+1)*vol:
						#raise Exception('ERROR: the condition \'cut<=const*log(n+1)*vol\' is not achieved for any r<1/c---lhs: %s, rhs: %s' %(cut,const*log(n+1)*vol))
						print('WARNING: the condition \'cut<=const*log(n+1)*vol\' is not achieved for any r<1/c---lhs: %s, rhs: %s' %(cut,const*log(n+1)*vol))
					break
				"""
				#######################
				#######################
				#######################
				## DEBUG:
				incr_cut = _incremental_cut(ball,new_ball_vertices,remaining_vertices,graph)
				cut_check = _cut(ball,remaining_vertices,graph)
				#cut = cut_check
				tolerance = 0.0000000001
				if abs(cut-cut_check) > tolerance:
					print('ERROR! cut incremental: %s, cut from scratch: %s' %(str(cut),str(cut_check)))
					sys.exit()
				#######################
				#######################
				#######################
				"""
			clusters.append(ball)
			for v in ball:
				remaining_vertices.remove(v)
	return clusters

#rounding algorithm proposed in Charikar et al., "Clustering with Qualitative Information", JCSS 2005
def round_charikar(x,id2vertexpair,id2vertex,edges,graph,lp_cost):
	clusters = []
	n = len(id2vertex)
	remaining_vertices = set(id2vertex.keys())
	shuffled_pairs = [id2vertexpair[h] for h in range(len(x)) if x[h]>2/3]
	random.shuffle(shuffled_pairs)
	#r = 1/3 - eps

	for (i,j) in shuffled_pairs:
		if i in remaining_vertices and j in remaining_vertices:
			ball = {i}
			r = _best_radius_charikar(i,x,remaining_vertices,n,graph,lp_cost)
			for k in remaining_vertices:
				if k!=i and x[_vertex_pair_id(i,k,n)]<=r:
					ball.add(k)
			clusters.append(ball)
			for k in ball:
				remaining_vertices.remove(k)
	for i in remaining_vertices:
		clusters.append({i})

	return clusters

def _best_radius_charikar(ball_center,x,valid_vertices,n,graph,lp_cost):
	d = _sorted_distances(ball_center,valid_vertices,n,x)
	r = 1/3 - eps

	i = 0
	while i<len(d) and d[i][1]<1/3:
		i += 1
	d = d[0:i]

	if d:
		r = 0
		ball = {ball_center}
		cut = _cut(ball,valid_vertices,graph)
		vol = lp_cost/n + _vol(ball_center,ball,valid_vertices,graph,n,x,r)
		best_ratio = float('inf') if vol==0 else cut/vol
		while d:
			new_r = d[0][1] #d[0][1]: minimum distance in d
			i = 0
			while i<len(d) and d[i][1]<=new_r:
				i += 1
			new_ball_vertices = {v for (v,d) in d[0:i]}
			cut += _incremental_cut(ball,new_ball_vertices,valid_vertices,graph) #cut can be computed incrementally
			ball.update(new_ball_vertices)
			vol = lp_cost/n + _vol(ball_center,ball,valid_vertices,graph,n,x,new_r)
			ratio = float('inf') if vol==0 else cut/vol
			if ratio<best_ratio:
				best_ratio = ratio
				r = new_r
			d = d[i:]
	return r

#well-established randomized, linear-time algorithm for correlation clustering, achiving constant-factor approximation guarantees on complete graphs
#see Ailon et al., "Aggregating inconsistent informa- tion: Ranking and clustering", JACM 2008
def kwikcluster(id2vertex,graph):
	clusters = []
	n = len(id2vertex)
	remaining_vertices = set(id2vertex.keys())
	shuffled_vertices = list(id2vertex.keys())
	random.shuffle(shuffled_vertices)

	for u in shuffled_vertices:
		if u in remaining_vertices:
			cluster = {u}
			if u in graph:
				for v in graph[u]:
					if v in remaining_vertices:
						(wp,wn) = graph[u][v]
						if wp > wn:
							cluster.add(v)
			clusters.append(cluster)
			for v in cluster:
				remaining_vertices.remove(v)

	return clusters

def _CC_cost(clustering,graph):
	cost = 0
	vertex2cluster = {}
	cid = 0
	for cluster in clustering:
		for u in cluster:
			vertex2cluster[u] = cid
		cid += 1
	for u in graph.keys():
		for v in graph[u]:
			if u<v:
				(wp,wn) = graph[u][v]
				if vertex2cluster[u] == vertex2cluster[v]:
					cost += wn
				else:
					cost += wp
	return cost

def _lp_solution_cost(lp_var_assignment,graph,num_vertices):
	cost = 0
	for u in graph.keys():
		for v in graph[u]:
			if u<v:
				xuv = lp_var_assignment[_vertex_pair_id(u,v,num_vertices)]
				(wp,wn) = graph[u][v]
				if wp>wn:
					cost += wp*xuv
				else:
					cost += wn*(1-xuv)
	return cost

def _all_edgeweights_sum(graph):
	sum = 0
	for u in graph.keys():
		for v in graph[u]:
			if u<v:
				(wp,wn) = graph[u][v]
				sum += wp
				sum += wn
	return sum

def _all_negativeedgeweight_sum(graph):
	sum = 0
	for u in graph.keys():
		for v in graph[u]:
			if u<v:
				(wp,wn) = graph[u][v]
				if wn>wp:
					sum += wn
	return sum

def _max_edgeweight_gap(graph):
	maxgap = 0
	for u in graph.keys():
		for v in graph[u]:
			if u<v:
				(wp,wn) = graph[u][v]
				maxgap = max(maxgap,abs(wp-wn))
	return maxgap

#for DEBUG
def _check_constraint_correspondence(A,c_dict):
	if len(A) != len(c_dict):
		return False

	for i in range(len(A)):
		v = A[i]
		posA = [j for j in range(len(v)) if v[j] == 1]
		if len(posA) != 1:
			raise Exception('Malformed constraint on \'A\'')
		negA = [j for j in range(len(v)) if v[j] == -1]
		if len(negA) != 2:
			raise Exception('Malformed constraint on \'A\'')

		c = c_dict[i]
		cs = str(c)[:-5].replace('- ','-').replace('+ ','')
		monomials = cs.split()
		#print(monomials)
		#print(cs)
		posc = [int(s.split('_')[1]) for s in monomials if s[0] != '-']
		if len(posc) != 1:
			print(i)
			print('Malformed constraint on \'c\'---c: ' + str(c))
			raise Exception('Malformed constraint on \'c\'---c: ' + str(c))
		negc = [int(s.split('_')[1]) for s in monomials if s[0] == '-']
		if len(negc) != 2:
			print(i)
			print(negc)
			print('Malformed constraint on \'c\'---c: ' + str(c))
			raise Exception('Malformed constraint on \'c\'---c: ' + str(c))

		if posA[0] != posc[0] or set(negA) != set(negc):
			return False

	return True

#for DEBUG
def _check_objective_correspondence(c,objective):
	objective_tokens = str(objective).replace('- ','-').replace('+ ','+').split()
	monomials = {}
	for s in objective_tokens:
		s_tokens = s.split('*')
		coeff = float(s_tokens[0])
		var_index = int(s_tokens[1].split('_')[1])
		monomials[var_index] = coeff

	for i in range(len(c)):
		if c[i] == 0 and i in monomials.keys() and monomials[i] != 0:
			return False
		if c[i] != 0 and (i not in monomials.keys() or monomials[i] != c[i]):
			return False

	return True

#for DEBUG
def _check_clustering(clustering,num_vertices):
	vertex2cluster = {}
	cid = 0
	for cluster in clustering:
		for u in cluster:
			if u not in vertex2cluster:
				vertex2cluster[u] = set()
			vertex2cluster[u].add(cid)
		cid += 1

	for u in range(num_vertices):
		if u not in vertex2cluster or len(vertex2cluster[u]) != 1:
			return False

	return True