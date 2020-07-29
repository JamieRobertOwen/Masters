using CSV, DataFrames
#CSV for file writing
#DataFrames for DataFrame type


"""generates test data"""
#default variables
NumPatients = 20
numSpecialists = 6
appointmentLengths = [1,2,3]
numAppoints = 80
function junk()

    files = ["a.txt", "b.txt", "c.txt"]
    fvars = Dict()
    for (n, f) in enumerate(files)
        fvars["x_$(n)"] = f
    end


    #generate TestInfo file
    PatientInfo = [0; NumPatients]
    SpecialistInfo = [0; numSpecialists]
    AppointmentInfo = [1; appointmentLengths]
    TestCasesInfo = [0; numAppoints]

    TestInfo = DataFrame(PatientInfo = PatientInfo, SpecialistInfo = SpecialistInfo, AppointmentInfo = AppointmentInfo, TestCasesInfo = TestCasesInfo)
    CSV.write("testInfo.csv", TestInfo)
end


#generating test case
SpecialistID= rand(1:1:numSpecialists,numAppoints)
AppointmentLength = rand(appointmentLengths,numAppoints)
#only allowing patients to have one appointment
#PatientID=rand(1:1:numPatients,numAppoints)
PatientID= 1:1:numAppoints

#writing test case to file
Appointment = DataFrame(SpecialistID = SpecialistID,
    AppointmentLength = AppointmentLength,
    PatientID = PatientID)
CSV.write("testData.csv",Appointment)
