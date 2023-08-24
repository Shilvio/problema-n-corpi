timesBody = [[], [], [], [], [], []]
timeIter=[]
for i in range(6):

    #00+j+(i*10)
    
    times = [[], [], [], [], [], []]
    for j in range(5):
        times[0].append(1)
        times[1].append(2)
        times[2].append(3)
        times[3].append(4)
        times[4].append(5)
        times[5].append(6)
        if i==1 and j==4:
            times[0][2]=20
            times[1][2]=20
            times[2][2]=20
            times[3][2]=20
            times[4][2]=20
            times[5][2]=20
    timeIter.append(times)

#print(timeIter)

for bodys in timeIter:
    count=0
    for time in bodys:
        temp=0
        for exec in time:
            temp+=exec
        timesBody[count].append(temp/len(time))
        count+=1

print(timesBody)