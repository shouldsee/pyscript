def convert_to_intervals(t_obs,s_obs,interval,t_start,t_final):
    T=t_start
    t_new_obs=[]
    s_new_obs=[]
    i=0
    ti=t_obs[i]
    si=s_obs[i]
    while T<t_final:
        while ti<T and i<len(t_obs)-1:
            i+=1
            ti=t_obs[i]
            si=s_obs[i]
        t_new_obs.append(T)
        s_new_obs.append(si)
        T=T+interval
    return (t_new_obs,s_new_obs)
