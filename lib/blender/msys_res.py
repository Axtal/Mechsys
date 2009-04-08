import Blender
import msys_dict as di

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
        for idx, lbl in obj.properties['res'][s]['idx2lbl'].iteritems():
            res       = [r for r in obj.properties['res'][s][lbl]]
            mi, ma    = min(res), max(res)
            stat[lbl] = [mi, 'n%d'%res.index(mi), ma, 'n%d'%res.index(ma)]

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
               '  %6s:[min at (node#/elem#),  max at (node#/elem#)]'%'key',
               '  -----------------------------------------------------------------']
        # nodes
        for k, v in stat.iteritems():
            msg.append('  %6s:[min = %g at %s,  max = %g at %s]'%(k,v[0],v[1],v[2],v[3]))

        # extra
        if len(ext)>0:
            msg.append('\n')
            msg.append('  %6s:[max of max at (elem#)]'%'key')
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

    # extra nodes
    if obj.properties.has_key('res_nodes'):
        if not obj.properties['res_nodes']=='':
            arr = obj.properties['res_nodes'].split(',')
            nds = [int(n) for n in arr]
            for n in nds:
                f = open ('%s_nod%d.cal'%(obj.name,n),'w')
                # header
                lin = '%4s' % 'stage'
                for idx, lbl in obj.properties['res']['1']['idx2lbl'].iteritems():
                    lin = '%s  %4s' % (lin,lbl)
                lin += '\n'
                f.write(lin)
                # values
                for s, v in obj.properties['res'].iteritems():
                    lin = ' %4d  ' % int(s)
                    for idx, lbl in obj.properties['res'][s]['idx2lbl'].iteritems():
                        res = obj.properties['res'][s][lbl][n]
                        lin = '%s  %8.3e' % (lin,res)
                    lin += '\n'
                    f.write (lin)
                f.close()

    # stats
    f = open (obj.name+'.res','w')
    for s, v in obj.properties['res'].iteritems():
        msg = stage_stats (int(s), False, False)
        for line in msg: f.write (line+'\n')
        f.write ('\n')
    f.close()
