using JuMP, CSV, DataFrames, GLPK
#JuMP to construct problem
#GLPK as solver
#CSV for data imput and output
#DataFrames for type
""" first pass at implimenting MIP """

NumPatients = 20
numSpecialists = 6
appointmentLengths = [1,2,3]
numAppoints = 80
lengthOfDay = 50

TestInfo = CSV.read("testData.csv")


clinicMIP = Model(GLPK.Optimizer)


@variable(clinicMIP, strictlyLater[1:numAppoints,1:numAppoints], Bin)
@variable(clinicMIP, 0<=startingTime[1:numAppoints]<=lengthOfDay, Int)


@constraint(clinicMIP, constraint[j=1:2], sum(A[j,i]*x[i] for i=1:3) <= b[j])


#MVP - maximise patients scheduled

#constraints
#are two appointments simultaneous? and DVs
#simultaneous appointments - not same doctor
#simultaneous appointments - not same room
#patient scheduled at most once (in the MVP)
#appointment blocks must be consecutive
