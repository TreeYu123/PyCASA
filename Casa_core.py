#coding=utf-8
import xarray as xr
import os
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
os.getcwd()

'''
输入：潜在蒸散发pet (mm) ; 降水pcp (mm)
'''
def cal_netrad(pet,pcp):
    nrad=(pet*pcp)**0.5*(0.369+0.598*(pet/pcp)**0.5)
    return nrad
'''
输入：降水pcp (mm) ; 净辐射nrad(mj/m^2)
ref: 周广胜,张新时.自然植被净第一性生产力模型初探.植物生态学报,1995,03.
'''
def cal_et(pcp,nrad):
    aa=pcp*nrad*(pcp**2+nrad**2+pcp*nrad)
    bb=(pcp+nrad)*(pcp**2+nrad**2)
    et=aa/bb
    return et

'''计算水分胁迫
输入：潜在蒸散发pet (mm) ; 蒸发et (mm)
'''
def cal_wen(pet,et):
    wen=0.5+0.5*et/pet
    return wen

'''计算光能利用率
输入：一年的气温t2_ls (摄氏度) ; 一年的ndvi ndvi_ls; 一年的水分胁迫wen_ls; 最大光能效率 max_ue
最大光能效率ref: 朱文泉,潘耀忠,张锦水.中国陆地植被净初级生产力遥感估算[J].植物生态学报,2007(03):413-424.
'''
def cal_ten(t2_ls,ndvi_ls,wen_ls,max_ue):
    flag=[]
    for i in range(12):
        flag.append(float(ndvi_ls[i].mean().values))
    maxid=flag.index(max(flag))
    t_opt=t2_ls[maxid]   #植物生长的最适温度
    
    #t_e1=[];t_e2=[];
    ue=[]
    for i in range(12):
        ##计算在低温和高温时植物内在的生化作用对光合的限制而降低净第一性生产力
        met2=t2_ls[i].mean().values
        if met2<=-10:
            t_e1i=t_opt*0
        else:
            t_e1i=0.8+0.2*t_opt-0.005*t_opt**2
        #t_e1.append(t_e1i)
        ##计算环境温度从最适温度t_opt向高温或低温变化时植物光能利用率变小的趋势 t_e2
        if ((met2-t_opt).mean().values>10) or ((met2-t_opt).mean().values<-13):
            aa=1.184/(1+np.exp(0.2*(-10)))
            bb=1/(1+np.exp(0.3*(-10)))
            t_e2i=(aa*bb)/2*(t2_ls[i]*0+1)
        else:
            aa=1.184/(1+np.exp(0.2*(t_opt-10-t2_ls[i])))
            bb=1/(1+np.exp(0.3*(-t_opt-10+t2_ls[i])))
            t_e2i=aa*bb
        #t_e2.append(t_e2i)
        uei=t_e1i*t_e2i*wen_ls[i]*max_ue
        ue.append(uei)
    ####
    return ue

'''计算光合有效辐射
输入：太阳总辐射srad ((MJ*m2)/month) ; 植被ndvi  ndvi (modis取值0-10000)
'''
def cal_apar(srad,ndvi):
    sr=(1+(ndvi*0.0001))/(1-(ndvi*0.0001))  ##植被指数
    sr_max=sr.quantile(0.95).values
    #print(sr_max)
    
    fapr = (sr-1.08)/(sr_max-1.05) #植被层对入射光合有效辐射的吸收比例
    fapr.values[fapr.values>0.95]=0.95  ##fapr 计算依据？参考文献？
    
    apar=srad*fapr*0.5
    return apar
    
'''计算植被净初级生产力NPP
输入：光合有效辐射apar ; 光能利用率ue
'''
def cal_npp(apar,ue):
    npp=apar*ue
    npp.values[npp.values<0]=0
    return npp

    
def main():
    pth=os.getcwd()+'\\Output'
    pth1=os.getcwd()+'\\Input\\evp'
    pth2=os.getcwd()+'\\Input\\max_ue'
    pth3=os.getcwd()+'\\Input\\ndvi'
    pth4=os.getcwd()+'\\Input\\prcp'
    pth5=os.getcwd()+'\\Input\\srad'
    pth6=os.getcwd()+'\\Input\\t2'
    yr=2001;yr_mue=2000
    mon= range(1,13)
    t2_ls=[];ndvi_ls=[];wen_ls=[];srad_ls=[]
    nm2=pth2+'\\'+str(yr_mue)+'.tif'
    ds_mue0=xr.open_rasterio(nm2)
    ds_mue=ds_mue0.where((ds_mue0>-9999)&(ds_mue0<9999))
    for mi in mon:
        nm1=pth1+'\\'+'Ep0'+str(yr)+str(mi).zfill(2)+'.tif'
        
        nm3=pth3+'\\'+'ndvi'+str(yr)+str(mi).zfill(2)+'.tif'
        nm4=pth4+'\\'+str(yr)+str(mi).zfill(2)+'.tif'
        nm5=pth5+'\\'+'sol'+str(yr)+str(mi).zfill(2)+'.tif'
        nm6=pth6+'\\'+str(yr)+str(mi).zfill(2)+'.tif'
        
        ###以基于lucc计算的ue坐标为基准进行匹配，避免出现矩阵大小不匹配而报错
        ds_evp0=xr.open_rasterio(nm1)
        ds_evp1=ds_evp0.interp(x=ds_mue.x.values,y=ds_mue.y.values,method='nearest')
        ds_evp=ds_evp1.where((ds_evp1>-9999)&(ds_evp1<9999))
        ds_ndvi0=xr.open_rasterio(nm3)
        ds_ndvi1=ds_ndvi0.interp(x=ds_mue.x.values,y=ds_mue.y.values,method='nearest')
        ds_ndvi=ds_ndvi1.where((ds_ndvi1>-9999)&(ds_ndvi1<9999))
        ds_prcp0=xr.open_rasterio(nm4)
        ds_prcp1=ds_prcp0.interp(x=ds_mue.x.values,y=ds_mue.y.values,method='nearest')
        ds_prcp=ds_prcp1.where((ds_prcp1>-9999)&(ds_prcp1<9999))
        ds_srad0=xr.open_rasterio(nm5)
        ds_srad1=ds_srad0.interp(x=ds_mue.x.values,y=ds_mue.y.values,method='nearest')
        ds_srad=ds_srad1.where((ds_srad1>-9999)&(ds_srad1<9999))
        ds_t20=xr.open_rasterio(nm6)
        ds_t21=ds_t20.interp(x=ds_mue.x.values,y=ds_mue.y.values,method='nearest')
        ds_t2=ds_t21.where((ds_t21>-9999)&(ds_t21<9999))
        #print(ds_t2)
        ######计算开始
        ds_nrad=cal_netrad(ds_evp,ds_prcp)
        ds_et=cal_et(ds_prcp,ds_nrad)
        ds_wen=cal_wen(ds_evp,ds_et)

        t2_ls.append(ds_t2)
        ndvi_ls.append(ds_ndvi)
        wen_ls.append(ds_wen)
        srad_ls.append(ds_srad)
        ##
    ###计算一年的npp g C*m^-2############################################
    ds_ue_ls=cal_ten(t2_ls,ndvi_ls,wen_ls,ds_mue)
    ds_npp=[];tm=[]
    for i,ds_ue in enumerate(ds_ue_ls):
        apar_i=cal_apar(srad_ls[i],ndvi_ls[i])
        ds_nppi=cal_npp(apar_i,ds_ue)
        ds_npp.append(ds_nppi)
        tm.append(dt.datetime(yr,i+1,1))
    
    ds_npp_all=xr.concat(ds_npp,dim='band')
    ds_npp_ulta=ds_npp_all.assign_coords(band=tm)
    ds_npp_ulta.name='npp'
    ds_npp_ulta.attrs={'unit':'g C*m^-2','long_name':'Net primary productivity'}
    outnm=pth+'\\npp'+str(yr)+'.nc'
    ds_npp_ulta.to_netcdf(outnm)
    print('Finish!')
    ##############################################################
    
    ###绘图测试1
    plt.rcParams['font.size'] = 8
    plt.rcParams['font.sans-serif']='times new roman'
    fig,ax=plt.subplots(1,1,figsize=(6.,4.1))
    lg=ds_npp_ulta.mean(axis=0).plot(ax=ax,add_colorbar=False)
    fig.colorbar(lg,extendrect=False,extend='both',label='NPP (g C*m$^{-2}$)',fraction=0.1)
    ax.set_title('annual mean %s'%yr)
    plt.savefig('annual_mean.png',dpi=600,bbox_inches='tight')

    ###绘图测试2
    fig,ax=plt.subplots(1,1,figsize=(6.,3.5))
    ds_npp_ulta.mean(axis=(1,2)).plot(ax=ax)
    ax.set_title('NPP mean %s'%yr)
    ax.set_xlabel('time')
    ax.set_ylabel('NPP (g C*m$^{-2}$)')
    plt.savefig('region_mean.png',dpi=600,bbox_inches='tight')
    #plt.show()
    #break
    

if __name__=='__main__':
    main()
