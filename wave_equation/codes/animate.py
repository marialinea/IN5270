from PIL import Image
import glob
import os

# Create the frames
frames = []
imgs = glob.glob("*.png")
sorted_imgs = []
for i in range(len(imgs)):
    string = imgs[i]
    string = string.split(".")
    sorted_imgs.append(string[0])
sorted_imgs = sorted(sorted_imgs)
sorted_imgs = [i + ".png" for i in sorted_imgs]

for i in sorted_imgs:
    new_frame = Image.open(i)
    frames.append(new_frame)

# Save into a GIF file that loops forever
frames[0].save('CosineHill.gif', format='GIF',
               append_images=frames[1:],
               save_all=True,
               duration=300, loop=0)
os.system("rm -rf *.png*")
