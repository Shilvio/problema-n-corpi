
def ceck(st,cm,who):
    st = st.split("\n")
    cm = cm.split("\n")
    for i in range(len(st)):
        if st[i]!=cm[i]:
            st2 = st[i].split(",")
            cm2 = cm[i].split(",")
            for j in range(len(st2)):
                if st2[j]!=cm2[j]:
                    st3 = st2[j]
                    cm3 = cm2[j]
                    for k in range(len(st3)):
                        if ((k != 8)and(k != 7) and (k != 6) and (k!=5)) and (st3[k] != cm3[k]):
                            print(who)
                            print(f"diverità alla riga {i+1}")
                            print(f"\tposizione {j+1}")
                            print(f"\t\t{st2[j]}\n\t\t{cm2[j]}")
                            break

def exec(s1,s2,c1,c2,m1,m2):
    stE = open("./St_Exaustive.txt","r").read()
    stB = open("./ST_Barnes-hut.txt","r").read()
    cudaE = open("./CUDA_Exaustive.txt","r").read()
    cudaB = open("./CUDA_Barnes-hut.txt","r").read()
    mpiE = open("./MPI_Exaustive.txt","r").read()
    mpiB = open("./MPI_Barnes-hut.txt","r").read()

    if stE!=cudaE:
        ceck(stE,cudaE,"CUDA exaustive")

    if stE!=mpiE:
        ceck(stE,mpiE,"MPI exaustive")

    if stB != cudaB:
        ceck(stB,cudaB,"CUDA Barnes-hut")

    if stB != mpiB:
        ceck(stB,mpiB,"MPI Barnes-hut")


exec(True,True,True,True,True,True)

