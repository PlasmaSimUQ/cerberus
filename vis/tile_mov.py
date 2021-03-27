# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 15:07:12 2016

@author: daryl
"""

import os, shutil
from multiprocessing import Pool
import subprocess, shlex
from PIL import Image, ImageFont, ImageDraw

#==============================================================================
#
#==============================================================================


def listdir_fullpath(d):
    return [os.path.join(d, f) for f in os.listdir(d)]

def tile(q):
    """
    tile images
    """

    out = q["out"]
    order = q["order"]

    qq = q["q"]

    nrows = 0
    ncols = 0
    for place in order:
        nrows = max(nrows, place[0])
        ncols = max(ncols, place[1])

    nrows += 1
    ncols += 1

    images = []
    widths = []
    heights = []

    for i, I in enumerate(q["set"]):
        
        if not I["image"]:
            images.append(None)
            continue

        img = Image.open(I["image"])
        ok = I["ok"]

        if "crop" in qq["dirs"][i]:
            [[x0, x1],[y0,y1]] = qq["dirs"][i]["crop"]
            img = img.crop((x0, y0, x1, y1))

        if "crop_relative" in qq["dirs"][i]:
            [[x0, x1],[y0,y1]] = qq["dirs"][i]["crop_relative"]
            w, h = img.size
            x0 = int(x0*w)
            x1 = int(x1*w)
            y0 = int(y0*h)
            y1 = int(y1*h)
            img = img.crop((x0, y0, x1, y1))

        if qq["max_individual_size"] != None:
            img.thumbnail(qq["max_individual_size"])

        if "label_color" in qq["dirs"][i]:
            label_color = qq["dirs"][i]["label_color"]
        else:
            label_color = (0,0,0)

        if "label" in  qq["dirs"][i]:
            label = qq["dirs"][i]["label"]
            draw = ImageDraw.Draw(img)
            font = ImageFont.truetype(qq["font"], qq["label_size"])
            text_size = font.getsize(label)
            
            w, h = img.size
            bnd = 100
            if "label_loc" in qq["dirs"][i]:
                loc = qq["dirs"][i]["label_loc"]
                if loc[0] == "+":
                    x = w - text_size[0]/2.0 - bnd
                elif loc[0] == "-":
                    x = text_size[0]/2.0 + bnd
                elif loc[0] == "o":
                    x = w/2.0
                elif type(loc[0]) != str:
                    x = loc[0]*w
                
                if loc[1] == "+":
                    y = text_size[1]/2.0 + bnd
                elif loc[1] == "-":
                    y = h - text_size[1]/2.0 - bnd
                elif loc[1] == "o":
                    y = h/2.0
                elif type(loc[1]) != str:
                    y = (1-loc[1])*h
            else:
                x = text_size[0]/2.0 + bnd
                y = text_size[1]/2.0 + bnd
            
            x = int(x)
            y = int(y)
                
            
            draw.text((x, y),label,label_color,font=font)

        if qq["tag_failed"]:
            if not ok:
                W, H = img.size
                draw = ImageDraw.Draw(img)
                font = ImageFont.truetype(qq["font"], qq["label_size"])
                txt = "FAILED"
                w, h = draw.textsize(txt)
                draw.text(((W-w)/2,(H-h)/2),txt,label_color,font=font)
                img = img.convert('1')

        images.append(img)
        w, h = img.size
        widths.append(w)
        heights.append(h)


    width = max(widths)
    height = max(heights)

    total_height = nrows*height
    total_width = ncols*width

    size = [total_width, total_height]

    # make sure we have frames that are divisible by 2
    for i, s in enumerate(size):
        ns = 2*int(s/2.0)
        if s != ns:
            size[i] += 1

    img = Image.new('RGB', size, color=(255,255,255,0))

    for i, im in enumerate(images):
        if not im:
            continue
        x_offset = order[i][1]*width
        y_offset = order[i][0]*height
        img.paste(im, (x_offset, y_offset))
        #im.close()

    img.save(out)
    #img.close()

    print(os.path.split(out)[1])

    return

def tile_movie(q):
    """
    given a list of folders join all the images in those folders together
    """


    items = {"only_frames":False, "font":"FreeMono.ttf"}

    for item, value in items.items():
        if item not in q:
            q[item] = value

    if "tag_failed" not in q:
        q["tag_failed"] = False

    q["img_out"] = os.path.join(q["out_dir"], q["name"])
    q["mov_out"] = q["out_dir"]

    # make sure the output directory exists
    try:
        os.makedirs(q["img_out"])
    except OSError:
        if not os.path.isdir(q["img_out"]):
            raise

    try:
        os.makedirs(q["mov_out"])
    except OSError:
        if not os.path.isdir(q["mov_out"]):
            raise

    # get a list of all images we want to merge
    order = []
    imgs = []
    levels = []
    clean = []
    for d in q["dirs"]:
        
        if not os.path.exists(d["dir"]):
            imgs.append([])
            order.append(d["tile"])
            continue
        
        # if we have passed in a movie then get the frames
        if not os.path.isdir(d["dir"]):
            d["dir"] = extract_frames(d["dir"], q)
            clean.append(d["dir"])
            
            
        print(d["dir"])

        levels.append(d["level"])
        sub_imgs = []
        for img in listdir_fullpath(d["dir"]):
            if ".png" not in img:
                continue
            num_str = img[-9:-4]
            number = int(num_str)

            if "trim" in q:
                if (number < q["trim"][0]) or (number > q["trim"][1]):
                    continue
            
            if "level" in os.path.split(img)[1]:
                tag = ["level-%i"%d["level"],
                       "level_%i"%d["level"],
                       "level=%i"%d["level"]]
    
                for t in tag:
                    if t in img:
                        sub_imgs.append(img)
                        break
            else:
                sub_imgs.append(img)
                
        sub_imgs = sorted(sub_imgs)

#        if not sub_imgs:
#            raise RuntimeError("No images found")

        imgs.append(sub_imgs)
        order.append(d["tile"])

    # update name to include level id
    if q["level_label"]:
        q["name"] += "-level=%s"%(str(levels))

    q["name"] = q["name"].replace(" ","")

    # number of tiles
    n = len(imgs)
    # get the maximum frames
    N = max([len(img_set) for img_set in imgs])

    last_good = []
    for j in range(n):
        last_good.append(0)

    # form them into sets
    image_sets = []
    for i in range(N):
        img_set = []
        for j in range(n):
            if not imgs[j]:
                img_set.append({"image":None, "ok":True})
                continue
            try:
                img_set.append({"image":imgs[j][i], "ok":True})
                last_good[j] = i
            except:
                img_set.append({"image":imgs[j][last_good[j]], "ok":False})

        d = {}
        d["set"] = img_set
        d["out"] = os.path.join(q["img_out"], q["name"]+"-%05d.png"%i)
        d["order"] = order
        d["q"] = q

        image_sets.append(d)


    # clean out the destination folder
    files = listdir_fullpath(q["img_out"])
    for f in files:
        if q["name"] in f:
            os.remove(f)

    # merge them
    if q["cores"] == 1:
        for iset in image_sets:
            tile(iset)
    else:
        pool = Pool(processes=q["cores"])
        pool.imap_unordered(tile, image_sets)
        pool.close()
        pool.join()

    out_name = os.path.join(q["mov_out"], q["name"]) + ".mp4"

    try:
        proc=subprocess.Popen(shlex.split("rm %s"%out_name))
    except:
        pass

#    call = "avconv -framerate %i -f image2 -i "%q["framerate"]
#    call += os.path.join(q["img_out"], q["name"]+"-%05d.png")
#    call += " "
#    call += out_name

    if not q["only_frames"]:
        prm = {"framerate":q["framerate"],
               "pattern":os.path.join(q["img_out"], q["name"]+"-%05d.png"),
               "output":out_name
               }
    
        call = "ffmpeg -framerate {framerate} -i {pattern} -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p {output}".format(**prm)
    
    
    
    
        proc=subprocess.Popen(shlex.split(call))
        proc.communicate()


    # clean up
    for f in clean:
        shutil.rmtree(f)


    return

def tile_all(Q):
    for q in Q:
        tile_movie(q)
    return

def extract_frames(mov, q):
    """
    mov: path to a movie file
    extract all the frames from a movie
    """

    path, name = os.path.split(mov)
    name, ext = os.path.splitext(name)

    new = os.path.join(path, name)

    # make sure the output directory exists
    try:
        os.makedirs(new)
    except OSError:
        if not os.path.isdir(new):
            raise


    # extract frames

    prm = {"mov":mov,
        "framerate":q["framerate"],
        "pattern":os.path.join(new, name+"-%05d.png"),
        }

    call = "avconv -i {mov} -q:v 1 -r {framerate} -f image2 {pattern}".format(**prm)


    proc=subprocess.Popen(shlex.split(call))
    proc.communicate()


    return new

#==============================================================================
#
#==============================================================================

if __name__ == "__main__":
    
    mov = "/home/uqdbond1/CODE/cerberus/Exec/movie-level=0.mpg"
    path = extract_frames(mov, {"framerate":10})


    q = {}
    q["dirs"] = []
    q["dirs"].append({"dir":path, "level":0, "tile":[0,0], "label":"TRMI-L1-D0.1"})

    q["name"] = "TRMI-L1-D0.1"

    q["out_dir"] = "/home/uqdbond1/CODE/cerberus/Exec/movie"

    q["trim"] = [0, 999999]
    q["crop"] = False #[[180,1150],[0,530]]
    q["max_individual_size"] = (1920, 1080)
    q["label_size"] = 16
    q["framerate"] = 30
    q["cores"] = 10
    q["level_label"] = False


    tile_movie(q)
    
    print("DONE")
