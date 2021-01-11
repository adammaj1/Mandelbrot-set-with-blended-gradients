How to use [blend modes](https://en.wikipedia.org/wiki/Blend_modes) to make better images ( of Mandelbrot sets and ...)  
These images are inspired by great technique used by Adam Saka in [KFMovieMaker](https://en.wikibooks.org/wiki/Fractals/kallesfraktaler#KFMovieMaker) and his videos 

# Blend modes

Blend types ( modes) 

```c

typedef enum  {
	no = 0,
	average = 1	
		} BlendType; 
```




How to compute blended color from 2 iput colors ? Here color is a single number (shades of gray)  

```c
unsigned char GiveBlendedColor( double c1, double c2, BlendType Blend){

	unsigned char t;
	
	switch (Blend){
	
		case average: {t = (c1+c2)/2.0; break;}
		
		
		default: {}
	
	}
	
	return  t;


}

```


# examples


## potential (Linear step) plus slope

Here exterior of the Mandelbrot set is:
* described by potential
* coloured with 1D gray gradients

Interior is solid colour blue


First image is slope ( normal map = [Lambert reflection (Illumination model or shader )](https://gitlab.com/adammajewski/mandelbrot_wiki_ACh#using-complex-potential-and-lambert-reflection-illumination-model-or-shader-)   
![](./images/normal.png)  


Second image: 
* level sets of potential
* [linear step function](https://github.com/adammaj1/1D-RGB-color-gradient#gray-linear-colormap)

![](./images/Linear_steps.png "Gray Linear gradient ( colormap)")  
![](./images/Linear_steps_2D.png "RGB profiles of the GrayL colormap")    

![](./images/step_linear.png)  


3-rd image = Result of blending above images in average mode  

![](./images/average.png)  



## Potential ( sqrt step ) + slope


First image is the same as above = ( normal map = [Lambert reflection (Illumination model or shader )](https://gitlab.com/adammajewski/mandelbrot_wiki_ACh#using-complex-potential-and-lambert-reflection-illumination-model-or-shader-)



Second image: 
* level sets of potential
* [sqrt step function](https://github.com/adammaj1/1D-RGB-color-gradient#gray-sqrt-colormap). It is inverterd

```c
 p = 1.0 - p;
```



![](./images/Sqrt_steps.png "Gray Linear gradient ( colormap)")  
![](./images/Sqrt_steps_2D.png "RGB profiles of the GrayL colormap")    

![](./images/step_sqrt.png)  


3-rd image = Result of blending above images in average mode  

![](./images/average_sqrt.png)  



# see also
* [KFMovieMaker](https://www.maths.town/after-effects-plugins/kfmoviemaker/kfmoviemaker-download-and-installation) by Adam Saka
  * [github repo](https://github.com/adamsaka/KFMovieMaker)
  * [in wikibooks](https://en.wikibooks.org/wiki/Fractals/kallesfraktaler#KFMovieMaker)
* [color-blend](https://github.com/loilo/color-blend) by Florian Reuschel



# program
c console program for creating static images (png)

## dependencies
* gcc
* bash
* Image Magic convert ( from ppm to png )


## to run

```
make
```





# licence
[LICENCE](LICENCE)


# git
```git
git init
git add README.md
git commit -m "first commit"
git branch -M main
git remote add origin git@github.com:adammaj1/Mandelbrot-set-with-blended-gradients.git
git push -u origin main
```

local repo:  ~/Mandelbrot-set-with-blended-gradients

