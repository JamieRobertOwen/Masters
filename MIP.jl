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
strictSlack = 1
maxOvertime = 5


#TestInfo = CSV.read("testData.csv")
TestInfo = DataFrame!(CSV.File("testData.csv"))


#clinicMIP = Model(with_optimizer(CPLEX.Optimizer))
clinicMIP = Model(CPLEX.Optimizer)

#set_optimizer_attribute(clinicMIP,"CPXPARAM_MIP_Display",3)
#5% optimality gap
set_optimizer_attribute(clinicMIP,"CPX_PARAM_EPGAP",0.1)
set_time_limit_sec(clinicMIP, 600)

magicA =2*cos(pi/8)/(1+cos(pi/8))
magicB =2*sin(pi/8)/(1+cos(pi/8))

SpecialistID=TestInfo[:,:SpecialistID]
AppointmentLength=TestInfo[:,:AppointmentLength]
AppointmentStd=TestInfo[:,:AppointmentStd]

#littleMstd=minimum(AppointmentStd)*0.9

maxAppointsInDay = floor(Integer,lengthOfDay/minimum(AppointmentLength))
bigMStd=maximum(AppointmentStd)*maxAppointsInDay
distanceApprox = zeros(maxAppointsInDay+1, maxAppointsInDay)


#for only 1 appointment return its std
aplusb = zeros(Integer,maxAppointsInDay,maxAppointsInDay)
b = zeros(Integer,maxAppointsInDay)

#work out b
rem=0
for i=2:maxAppointsInDay
    div = i
    count = 0
    while true
        rem = div % 2
        if rem == 0
            count += 1
        end
        div = div ÷ 2
        div = div + rem
        div == 1 && break
    end
    b[i]=count
end


#work out aplusb
for i =2:maxAppointsInDay
    maxVal = ceil(Integer,log(2,i))
    for j=1:(i-1)
        aplusb[i,j] = min(aplusb[i-1,j]+1,maxVal)
    end
    aplusb[i,i] = b[i]
end

#restructure b
b=hcat([vcat(zeros(Integer,i-1,1), ones(Integer,maxAppointsInDay-i+1)*b[i]) for i=1:maxAppointsInDay]...)


#work out distanceApprox
distanceApprox[2:(maxAppointsInDay+1),:]=(aplusb-b)*magicA +b*magicB
distanceApprox[2,1]=1

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

""" test variable """
#@variable(clinicMIP, 0<=stdRanked[1:numSpecialists,1:maxAppointsInDay])
@variable(clinicMIP, stdIndex[1:numSpecialists,1:maxAppointsInDay,1:numAppoints], Bin)
#@variable(clinicMIP, stdIsZero[1:numSpecialists,1:maxAppointsInDay],Bin)
@variable(clinicMIP, stdAppointNbr[1:numSpecialists,1:(maxAppointsInDay+1)], Bin)

@variable(clinicMIP,std[1:numSpecialists])

@objective(clinicMIP, Max, sum(beingScheduled))

""" test constraint """
function junk()
    #@constraint(clinicMIP, stdRanking[i=1:numSpecialists, j=1:maxAppointsInDay, k= 1:maxAppointsInDay; j<k],
    #stdRanked[i,j]>= stdRanked[i,k]) #when ranked lower indexes are higher

    #@constraint(clinicMIP, stdIndexTrue[i=1:numSpecialists,j=1:maxAppointsInDay,k=1:numAppoints],
    #stdRanked[i,j] >= stdIndex[i,j,k]*AppointmentStd[k]+(beingTreatedBy[k,i]-1)*bigMStd)
    #make sure that if scheduled stdRanked has to be greater than stdIndex

    #@constraint(clinicMIP, stdIndexEnsure[i=1:numSpecialists,j=1:maxAppointsInDay]
    #sum(stdIndex[i,j,k] for k =1:numAppoints) == numAppoints-j+1)
    #limit number of cheats stdIndex gives

    #@constraint(clinicMIP, stdIsZero[i=1:numSpecialists,j=1:maxAppointsInDay],
    #0<=(1-stdIsZero[i,j])*bigMstd-stdRanked[i,j])
    #need to change this so that when stdIsZero is 1 you get bigMStd if 0 you get littleMstd
end

@constraint(clinicMIP, stdRanking[i=1:numSpecialists, j=1:maxAppointsInDay, k= 1:maxAppointsInDay; j<k],
sum(AppointmentStd[l]*(stdIndex[i,j,l]-stdIndex[i,k,l]) for l=1:numAppoints)>= 0)
#when ranked lower indexes are higher

@constraint(clinicMIP, stdIndexEnsure[i=1:numSpecialists,j=1:maxAppointsInDay],
sum(stdIndex[i,j,k] for k=1:numAppoints) <= 1)
#each appointment can only have one index maximum

@constraint(clinicMIP, stdNumberOfZeros[i=1:numSpecialists],
sum(stdIndex[i,j,k] for j=1:maxAppointsInDay, k=1:numAppoints)
== sum(stdAppointNbr[i,j]*(j-1) for j=1:(maxAppointsInDay+1)))
#counts number of appoints specialist could have compared to maximum

@constraint(clinicMIP, stdIndexTrue[i=1:numSpecialists,k=1:numAppoints],
sum(stdIndex[i,j,k] for j=1:maxAppointsInDay)>=beingTreatedBy[k,i])
#if an appointment is scheduled must have an index

@constraint(clinicMIP, stdAppointNbrRepeats[i=1:numSpecialists],
sum(stdAppointNbr[i,j] for j=1:maxAppointsInDay) ==1)
#make sure that appointnumber is mapped correctly

@constraint(clinicMIP,stdApprox[i=1:numSpecialists, k=1:(maxAppointsInDay+1)],
sum(stdIndex[i,j,l]*distanceApprox[k,j]*AppointmentStd[l]
 for j=1:maxAppointsInDay, l=1:numAppoints)-stdAppointNbr[i,k]*bigMStd<=std[i])
#calculates the approximate standard deviation

@constraint(clinicMIP, overtime[i=1:numSpecialists],
sum(beingTreatedBy[j,i]*AppointmentLength[j] for j=1:numAppoints)+2*std[i] <= lengthOfDay+maxOvertime
)


@constraint(clinicMIP, ifScheduled[i=1:numAppoints],
 sum(beingTreatedBy[i,j] for j=1:numSpecialists) == beingScheduled[i])

@constraint(clinicMIP, Overlap1[i=1:numAppoints, j = 1:numAppoints; i != j],
 startingTime[i]-endTime[j]+lengthOfDay*(2-beingScheduled[i]-beingScheduled[j])-strictSlack>= -lengthOfDay*startsBefore[i,j])

@constraint(clinicMIP, Overlap2[i=1:numAppoints, j=1:numAppoints; i!=j],
endTime[j]-startingTime[i]>=lengthOfDay*(startsBefore[i,j]-1))

@constraint(clinicMIP, SurgeonCloning[i=1:numAppoints,j=1:numAppoints, k=1:numSpecialists; i!= j],
beingTreatedBy[i,k]+beingTreatedBy[j,k]<=3-startsBefore[i,j]-startsBefore[j,i])

@constraint(clinicMIP, runTime[i=1:numAppoints],
endTime[i]==startingTime[i]+AppointmentLength[i])

@constraint(clinicMIP, correctSpecialist[i=1:numAppoints,j=1:numSpecialists],
beingTreatedBy[i,j]<=specialistTable[i,j])


optimize!(clinicMIP)

startingTimeSol=value.(startingTime)
specialistArraySol = value.(beingTreatedBy)

#specialistSol = mapslices(x->findfirst(x),specialistArraySol.==1,dims=[1])

#specialistSol = vcat([findfirst(specialistArraySol[j,:].==1) for j=1:numAppoints]...)
specialistSol = vcat([findfirst(specialistArraySol[j,:].==1) for j=1:numAppoints]...)

#specialistSolTest=replace(specialistSol, nothing => 0)

#nothinginds = findall(isnothing,specialistSol)
#startingNothing = startingTimeSol[nothinginds
startingTimeSol =convert(Vector{Union{Nothing, Integer}},startingTimeSol)
startingTimeSol[findall(isnothing,specialistSol)].=nothing
#@constraint(clinicMIP, constraint[j=1:2], sum(A[j,i]*x[i] for i=1:3) <= b[j])


#MVP - maximise patients scheduled

#constraints
#are two appointments simultaneous? and DVs
#simultaneous appointments - not same doctor
#simultaneous appointments - not same room
#patient scheduled at most once (in the MVP)
#appointment blocks must be consecutive
