#using JuMP, CSV, DataFrames, GLPK
using JuMP, CSV, DataFrames, CPLEX
#JuMP to construct problem
#GLPK as solver
#CSV for data imput and output
#DataFrames for type
""" first pass at implimenting MIP """

#NumPatients = 80
numSpecialists = 6
#appointmentLengths = [1,2,3]
numAppoints = 80
lengthOfDay = 50
bigM = 1000
strictSlack = 0.5

TestInfo = CSV.read("testData.csv")


clinicMIP = Model(with_optimizer(CPLEX.Optimizer))
#set_parameter(clinicMIP,"CPXPARAM_MIP_Display",1)
SpecialistID=TestInfo[:,:SpecialistID]
AppointmentLength=TestInfo[:,:AppointmentLength]
#PatientID=TestInfo[:,:PatientID]
#specialistTable=zeros(numAppoints,numSpecialists)
function array2pos(j)
    fizz=zeros(1,numSpecialists)
    fizz[j]=1
    return fizz
end
specialistTable  = vcat([array2pos(SpecialistID[j]) for j=1:numAppoints]...)

@variable(clinicMIP, beingScheduled[1:numAppoints],Bin)
@variable(clinicMIP, startsBefore[1:numAppoints,1:numAppoints], Bin)
@variable(clinicMIP, 0<=startingTime[1:numAppoints]<=lengthOfDay, Int)
@variable(clinicMIP, beingTreatedBy[1:numAppoints,1:numSpecialists],Bin)
@variable(clinicMIP, 0<=endTime[1:numAppoints]<=lengthOfDay,Int)

@objective(clinicMIP, Max, sum(beingScheduled))

@constraint(clinicMIP, ifScheduled[i=1:numAppoints],
 sum(beingTreatedBy[i,j] for j=1:numSpecialists) == beingScheduled[i])

@constraint(clinicMIP, Overlap1[i=1:numAppoints, j = 1:numAppoints; i != j],
 startingTime[i]-endTime[j]+lengthOfDay*(2-beingScheduled[i]-beingScheduled[j])+strictSlack>= -lengthOfDay*startsBefore[i,j])

@constraint(clinicMIP, Overlap2[i=1:numAppoints, j=1:numAppoints; i!=j],
endTime[j]-startingTime[i]>=lengthOfDay*(startsBefore[i,j]-1))

@constraint(clinicMIP, SurgeonCloning[i=1:numAppoints,j=1:numAppoints, k=1:numSpecialists; i!= j],
beingTreatedBy[i,k]+beingTreatedBy[j,k]<=3-startsBefore[i,j]-startsBefore[j,i])

@constraint(clinicMIP, runTime[i=1:numAppoints],
endTime[i]==startingTime[i]+AppointmentLength[i])

@constraint(clinicMIP, correctSpecialist[i=1:numAppoints,j=1:numSpecialists],
beingTreatedBy[i,j]<=specialistTable[i,j])

#optimize!(clinicMIP)


#@constraint(clinicMIP, constraint[j=1:2], sum(A[j,i]*x[i] for i=1:3) <= b[j])


#MVP - maximise patients scheduled

#constraints
#are two appointments simultaneous? and DVs
#simultaneous appointments - not same doctor
#simultaneous appointments - not same room
#patient scheduled at most once (in the MVP)
#appointment blocks must be consecutive
