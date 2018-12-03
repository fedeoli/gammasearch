# cta-gamma-ray-analysis
* Python 3 is required
* Install all the libraries listed in `requirements.txt` (`pip install -r requirements.txt`)

## Detect sources on a single image

1. Run `python detect.py filepath --config configpath` from the cta-gamma-ray-analysis **root folder**
```
usage: detect.py [-h] [--config CONFIG] filepath

Detect Sources

positional arguments:
  filepath         Path of the fits file

optional arguments:
  -h, --help       show this help message and exit
  --config CONFIG  Path of config file (default = data/default.conf)
 ```
## Parameter tuning

1. Run `run.sh` from the cta-gamma-ray-analysis **root folder**
2. **Keep focus on displayed images**
3. You can change:
	- **Adaptive Threshold**
		- **blockSize**, Size of a pixel neighbourhood that is used to calculate a threshold value for the pixel
		- **const**, Constant subtracted from the mean or weighted mean
	- **Filters**
		- **number of median filter iterations**
		- **median filter kernel size**
		- **number of gaussian filter iterations**
		- **gaussian filter kernel size**
	- **Local Linear Stretching**
		- **stretch kernel size**
		- **stretch step size**: step size at every sliding window iteration
		- **stretch min bins**: minimum value within which the current window is normalized
		
		
**Legenda**:

Type the following keys to change mode (remember: **keep focus on the displayed images, don't write anything on the terminal**):

- `t`	adaptive threshold
- `f`	gaussian & median filter
- `l`	local stretch

then type `[1 .. 4]` to select the parameter and use `a` and `d` to change it.

Then you can type:

- `r`	to print results 
- `v`	to print current values
- `p`	to visualize intermediate steps



- `w` and `s` to change map

- `enter` to print legenda
- `esc`	to quit

	
