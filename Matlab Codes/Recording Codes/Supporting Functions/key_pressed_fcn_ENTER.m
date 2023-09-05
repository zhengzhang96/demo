function key_pressed_fcn_ENTER(fig_obj,~)

key = get(fig_obj, 'CurrentKey');

if strcmp(key,'return')
    global keep %#ok<GVMIS>
    keep=0;
end

end