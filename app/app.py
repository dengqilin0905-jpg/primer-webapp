import streamlit as st
from auth_ui import login_box, get_current_user, logout
from profile_ui import profile_page
from backbones_ui import backbones_page
from genes_ui import genes_page
from design_ui import design_page
from history_ui import history_page
from templates_ui import templates_page
from config import load_config
from file_safety import ensure_temp_dir, cleanup_old_temp_files

st.set_page_config(page_title="Primer WebApp", layout="wide")
CFG = load_config()

# Optional site-wide access password (shared with invited users)
SITE_PASSWORD = CFG.app_password
if SITE_PASSWORD:
    if not st.session_state.get("site_authed", False):
        st.title("Access Required")
        st.write("This app is currently restricted. Please enter the access password.")
        pwd = st.text_input("Access password", type="password")
        if st.button("Enter"):
            if pwd == SITE_PASSWORD:
                st.session_state["site_authed"] = True
                st.rerun()
            else:
                st.error("Incorrect password.")
        st.stop()

if CFG.is_server:
    try:
        td = ensure_temp_dir(CFG.temp_dir)
        cleanup_old_temp_files(td, CFG.temp_file_retention_hours)
    except Exception:
        pass

login_box()
user = get_current_user()

if not user:
    st.title("Primer WebApp")
    st.info("Please login or register in the sidebar.")
    st.stop()

st.sidebar.write(f"Logged in as: **{user['username']}**")
if st.sidebar.button("Logout"):
    logout()
    st.rerun()

tabs = st.tabs(["Profile", "Backbones", "Genes (Batch Import)", "Templates", "Design (Batch)", "History"])
with tabs[0]:
    profile_page(user)
with tabs[1]:
    backbones_page(user, CFG)
with tabs[2]:
    genes_page(user, CFG)
with tabs[3]:
    templates_page(user, CFG)
with tabs[4]:
    design_page(user, CFG)
with tabs[5]:
    history_page(user, CFG)
