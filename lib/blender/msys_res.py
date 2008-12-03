import Blender
import msys_dict as di

def at_node(stage_id):
    edm, obj, msh = di.get_msh()
    verts = msh.verts.selected()
    if not obj.properties.has_key('res'):
        if edm: Blender.Window.EditMode(1)
        raise Exception('Please, run analysis first')
    if len(verts)==0:
        if edm: Blender.Window.EditMode(1)
        raise Exception('Please, select at least 1 vertex')
    if len(verts)>3:
        if edm: Blender.Window.EditMode(1)
        raise Exception('Please, select at most 3 vertices')
    s   = str(stage_id)
    msg = ['=== Stage # '+s+' =====']
    for vidx in verts:
        msg.append('--- Node # %d -----'%(vidx))
        for k, key in di.key('dfv').iteritems():
            if obj.properties['res'][s].has_key(key):
                val = obj.properties['res'][s][key][vidx]
                msg.append('  %s = %g'%(key,val))
    print 'Results:'
    for item in msg: print item
    Blender.Draw.PupBlock('Results:',msg)
    if edm: Blender.Window.EditMode(1)


def stage_stats(stage_id,console=True,popup=False):
    # get obj and msh
    edm, obj, msh = di.get_msh()

    # check
    if not obj.properties.has_key('res'):
        if edm: Blender.Window.EditMode(1)
        raise Exception('Please, run analysis first')

    # stags
    s = str(stage_id)
    if obj.properties['res'].has_key(s):

        # values at nodes
        stat = {}
        for k, key in di.key('dfv').iteritems():
            res       = [r for r in obj.properties['res'][s][key]]
            mi, ma    = min(res), max(res)
            stat[key] = [mi, 'n%d'%res.index(mi), ma, 'n%d'%res.index(ma)]

        # extra values
        ext = {}
        if obj.properties['res'][s].has_key('extra'):
            exts = ['N', 'M', 'V']
            for e in exts:
                key = 'max_'+e
                if obj.properties['res'][s].has_key(key):
                    v = obj.properties['res'][s][key]
                    ext[key] = [v[0], v[1]]

        # caption
        msg = ['===== Stage # '+s+' =========================================================','\n'
               '  %6s:[min = ??? at (node#/elem#),  max = ??? at (node#/elem#)]'%'key',
               '  -----------------------------------------------------------------']
        # nodes
        for k, v in stat.iteritems():
            msg.append('  %6s:[min = %g at %s,  max = %g at %s]'%(k,v[0],v[1],v[2],v[3]))

        # extra
        if len(ext)>0:
            msg.append('\n')
            msg.append('  %6s:[max of max = ??? at (elem#)]'%'key')
            msg.append('  --------------------------------------------')
        for k, v in ext.iteritems():
            msg.append('  %6s:[max of max = %g at e%d]'%(k,v[0],v[1]))

        # pop up
        if popup: Blender.Draw.PupBlock('Statistics:',msg)

        # output to console
        if console:
            print '\nStatistics:'
            for m in msg: print m
        else:
            if edm: Blender.Window.EditMode(1)
            return msg

    if edm: Blender.Window.EditMode(1)


def report():
    # get object
    obj = di.get_obj()

    # check
    if not obj.properties.has_key('res'):
        if edm: Blender.Window.EditMode(1)
        raise Exception('Please, run analysis first')

    # stats
    f = open (obj.name+'.res','w')
    for s, v in obj.properties['res'].iteritems():
        msg = stage_stats (int(s), False, False)
        for line in msg: f.write (line+'\n')
        f.write ('\n\n')

    # extra nodes
    if obj.properties.has_key('res_nodes'):
        arr = obj.properties['res_nodes'].split(',')
        nds = [int(n) for n in arr]
        lin = '\n %8s' % 'Node #'
        for k, key in di.key('dfv').iteritems(): lin = '%s  %12s' % (lin,key)
        lin += '\n'
        f.write (lin)
        for n in nds:
            res = obj.properties['res'][s][key][n]
            lin = ' %8d  ' % n
            for k, key in di.key('dfv').iteritems(): lin = '%s  %8.3e' % (lin,res)
            lin += '\n'
            f.write (lin)

    f.close()
