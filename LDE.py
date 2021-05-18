import math

def FindCeiling(solution, vector, minlen):
	i = 1
	while Length(SumVectors(solution, MultiplyVectorByNumber(vector, i))) < minlen:
		print(str(Length(SumVectors(solution, MultiplyVectorByNumber(vector, i)))) + " " + str(minlen))
		i += 1
	return i

def FindFloor(solution, vector, minlen):
	i = -1
	while Length(SumVectors(solution, MultiplyVectorByNumber(vector, i))) < minlen:
		i -= 1
	return i

def MinSolution(a, b):
	vectors = CreateSolutionMatrix(a, b)
	i = 0
	while vectors[i][0] == 0:
		i += 1
	gcd = vectors[i][0]
	solution = vectors[i][1 : len(vectors[i])]
	solution = MultiplyVectorByNumber(solution, b // gcd)
	minlen = Length(solution)
	parameters = []
	for j in range(len(vectors)):
		if (j != i):
			parameters.append(vectors[j][1 : len(vectors[j])])
	found, solution = BruteForce(a, b, solution, parameters, Length(solution))
	while (found):
		found, solution = BruteForce(a, b, solution, parameters, Length(solution))
	return solution

def Increase(parametersMultipliers, floors, ceilings):
	parametersMultipliers[0] += 1
	for i in range(len(parametersMultipliers)):
		if (parametersMultipliers[i] > ceilings[i]):
			if (i == len(parametersMultipliers) - 1):
				return True
			parametersMultipliers[i + 1] += 1
			parametersMultipliers[i] = floors[i]
	return False

def TotalSolution(solution,parameters,parametersMultipliers):
	totalSolution = solution
	for i in range(len(parameters)):
		totalSolution = SumVectors(totalSolution, MultiplyVectorByNumber(parameters[i],parametersMultipliers[i]))
	return totalSolution

def SumOfVector(vector):
	count = 0;
	for i in range(len(vector)):
		count += vector[i]
	return count

def BruteForce(a, b, solution, parameters, minlen):
	floors = []
	ceilings = []
	parametersMultipliers = []
	for j in range(len(parameters)):
		floors.append(FindFloor(solution, parameters[j], minlen))
		parametersMultipliers.append(floors[j])
		ceilings.append(FindCeiling(solution, parameters[j], minlen))
	ended = False
	while(not ended and not (Length(TotalSolution(solution, parameters, parametersMultipliers)) < minlen)):
		if (Increase(parametersMultipliers, floors, ceilings)):
			ended = True
	if (Length(TotalSolution(solution, parameters, parametersMultipliers)) < minlen):
		return True, TotalSolution(solution, parameters, parametersMultipliers)
	return False, solution

def AllSolutions(a, b):
	vectors = CreateSolutionMatrix(a, b)
	i = 0
	while vectors[i][0] == 0:
		i += 1
	gcd = vectors[i][0]
	solution = vectors[i][1 : len(vectors[i])]
	solution = MultiplyVectorByNumber(solution, b // gcd)
	for j in range(len(solution)):
		solution[j] = str(solution[j])
	for j in range(len(vectors)):
		if (vectors[j][0] == 0):
			for k in range(1, len(vectors[j])):
				solution[k - 1] += " + " + str(vectors[j][k]) + "*t" + str(j)
	return solution

def PartialSolution(a, b):
	vectors = CreateSolutionMatrix(a, b)
	i = 0
	while vectors[i][0] == 0:
		i += 1
	gcd = vectors[i][0]
	solution = vectors[i][1 : len(vectors[i])]
	solution = MultiplyVectorByNumber(solution, b // gcd)
	return solution
	

def CreateSolutionMatrix(a, b):
	vectors = []
	for i in range(len(a)):
		temp_vector = []
		temp_vector.append(a[i])
		for j in range(len(a)):
			if (i == j):
				temp_vector.append(1)
			else:
				temp_vector.append(0)
		vectors.append(temp_vector)
	for i in range(len(vectors)):
		if (vectors[i][0] < 0):
			vectors[i][0] = MultiplyVectorByNumber(vectors[i][0], -1)
	while not CheckVectors(vectors):
		minVector = FindMinVector(vectors)
		if (minVector is None):
			return None
		DivideVectorsByVector(vectors, minVector)
	return vectors

def CheckVectors(vectors):
	count = 0
	for i in range(len(vectors)):
		if (vectors[i][0] == 0):
			count += 1
	if (count == len(vectors) - 1):
		return True
	else:
		return False

def FindFirstNonNegativeVector(vectors):
	for i in range(len(vectors)):
		if (vectors[i] != 0):
			return vectors[i]
	return None

def FindMinVector(vectors):
	minVector = FindFirstNonNegativeVector(vectors)
	if (minVector is None):
		return None
	_min = minVector[0]
	for i in range(len(vectors)):
		if (vectors[i][0] < _min and vectors[i][0] != 0):
			_min = vectors[i][0]
			minVector = vectors[i]
	return minVector

def MultiplyVectorByNumber(vector, number):
	vec = []
	for i in range(len(vector)):
		vec.append(vector[i] * number)
	return vec

def SumVectors(vector1, vector2):
	vec = []
	for i in range(len(vector1)):
		vec.append(vector1[i] + vector2[i])
	return vec

def SubstractVectors(vector1, vector2):
	vec = []
	for i in range(len(vector1)):
		vec.append(vector1[i] - vector2[i])
	return vec

def DivideVectorsByVector(vectors, vector):
	for i in range(len(vectors)):
		if (vectors[i][0] != vector[0]):
			for j in range(len(vectors[i])):
				while vectors[i][0] >= vector[0]:
					vectors[i] = SubstractVectors(vectors[i], vector)

def gcd_extended(a, b):
	if a == 0:
		return (b, 0, 1)
	else:
		gcd, x, y = gcd_extended(b % a, a)
	return (gcd, y - (b // a) * x, x)

def GCD(a, b):
	if (b == 0):
		return a
	else:
		return GCD(b, a%b)

def GCD_multiple(a):
	if (len(a) == 2):
		return GCD(a[0],a[1])
	else:
		b = a[len(a) - 1]
		a.pop()
		return GCD(GCD_multiple(a), b)

def Length(a):
	count = 0
	for item in a:
		count += item**2
	return math.sqrt(count)

def CheckSolution(a,b,c):
	if (len(a) != len(b)):
		return False
	count = 0;
	for i in range(len(a)):
		count += a[i]*b[i]
	if (count == c):
		return True
	else:
		return False

def SolveLDEwith2vars(a, b, c):
	gcd = GCD(a, b)
	if (c % gcd != 0):
		return "Equation doesn't have integer solutions", 0
	else:
		a = a // gcd
		b = b // gcd
		c = c // gcd
		gcd2, x, y = gcd_extended(min(a,b), max(a,b))
		return x*c, y*c

def SolveAllLDEwith2vars(a, b, c):
	gcd = GCD(a, b)
	x0, y0 = SolveLDEwith2vars(a, b, c)
	if (x0 == "Equation doesn't have integer solutions"):
		return x0
	x = "x = " + str(x0) + " - n * " + str(b // gcd)
	y = "y = " + str(y0) + " + n * " + str(a // gcd)
	return x, y

def FindMinSolutionOfLDEWith2vars(a, b, c):
	gcd = GCD(a, b)
	x0, y0 = SolveLDEwith2vars(a, b, c)
	if (x0 == "Equation doesn't have integer solutions"):
		return x0
	x = x0
	y = y0
	n = 0
	while math.sqrt((x - (b // gcd) * n)**2 + (y + (a // gcd) * n)**2) < math.sqrt(x0**2 + y0**2):
		n += 1
	x = "x = " + str(x0 - (n - 1) * (b // gcd))
	y = "y = " + str(y0 + (n - 1) * (a // gcd))
	return x, y