from matplotlib.gridspec import GridSpec
import numpy as np
import tkinter as tk
from tkinter import messagebox, filedialog

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.patches import Circle


def compute_illuminance(
    W_mm, H_mm,
    Wres, Hres,
    xL_mm, yL_mm, zL_mm,
    I0,
    cx_mm, cy_mm, R_mm
):
    # координатная сетка в мм (для отображения)
    xs_mm = np.linspace(-W_mm / 2.0, W_mm / 2.0, Wres)
    ys_mm = np.linspace(-H_mm / 2.0, H_mm / 2.0, Hres)
    X_mm, Y_mm = np.meshgrid(xs_mm, ys_mm)

    # переводим всё в метры для физически корректного E (Вт/м^2)
    X = X_mm / 1000.0
    Y = Y_mm / 1000.0
    xL = xL_mm / 1000.0
    yL = yL_mm / 1000.0
    zL = zL_mm / 1000.0
    cx = cx_mm / 1000.0
    cy = cy_mm / 1000.0
    R = R_mm / 1000.0

    # расстояние до источника в метрах
    r2 = (X - xL) ** 2 + (Y - yL) ** 2 + zL ** 2

    # Ламбертовский источник:
    # I(θ) = I0 cosθ
    # E = I(θ) cosθ / r^2 = I0 cos^2θ / r^2
    # cosθ = zL / r => E = I0 * zL^2 / r^4
    # При r в метрах: E в Вт/м^2 (sr безразмерен)
    E = I0 * (zL ** 2) / (r2 ** 2)

    # маска круга в мм (геометрия на плоскости)
    circle_mask = (X_mm - cx_mm) ** 2 + (Y_mm - cy_mm) ** 2 <= R_mm ** 2
    circle_values = E[circle_mask]

    # нормировка только по кругу для визуализации
    E_norm = np.zeros_like(E, dtype=np.uint8)
    if circle_values.size > 0:
        Emax_c = float(circle_values.max())
        Emin_c = float(circle_values.min())
        Eavg_c = float(circle_values.mean())

        if Emax_c > 0:
            E_norm[circle_mask] = (
                circle_values / Emax_c * 255.0
            ).astype(np.uint8)
    else:
        Emin_c = Emax_c = Eavg_c = 0.0

    # контрольные точки в Вт/м^2
    def point_E(x_mm, y_mm):
        x = x_mm / 1000.0
        y = y_mm / 1000.0
        r2_pt = (x - xL) ** 2 + (y - yL) ** 2 + zL ** 2
        return float(I0 * (zL ** 2) / (r2_pt ** 2))

    control_points = {
        "Центр круга": point_E(cx_mm, cy_mm),
        "Пересечение с +X": point_E(cx_mm + R_mm, cy_mm),
        "Пересечение с -X": point_E(cx_mm - R_mm, cy_mm),
        "Пересечение с +Y": point_E(cx_mm, cy_mm + R_mm),
        "Пересечение с -Y": point_E(cx_mm, cy_mm - R_mm),
    }

    return {
        "E": E,                          # Вт/м^2
        "E_norm": E_norm,                # 0-255, только внутри круга
        "xs_mm": xs_mm,
        "ys_mm": ys_mm,
        "circle_mask": circle_mask,
        "circle_stats": {
            "E_min": Emin_c,
            "E_max": Emax_c,
            "E_avg": Eavg_c,
        },
        "control_points": control_points,
    }


# ----------------- GUI + встроенный matplotlib -----------------

root = tk.Tk()
root.title("ЛР3: Расчет освещенности от точечного источника (Ламбертовская диаграмма)")

# Основной фрейм: слева управление, справа графика
frame_controls = tk.Frame(root, padx=10, pady=10)
frame_controls.grid(row=0, column=0, sticky="n")

frame_plot = tk.Frame(root, padx=10, pady=10)
frame_plot.grid(row=0, column=1, sticky="nsew")

root.geometry("1400x900")
frame_plot.pack_propagate(False)

root.grid_columnconfigure(1, weight=1)
root.grid_rowconfigure(0, weight=1)

# ----- Поля ввода -----

# область
tk.Label(frame_controls, text="Ширина области W, мм").grid(row=0, column=0, sticky="w")
entry_W = tk.Entry(frame_controls, width=10)
entry_W.grid(row=0, column=1)
entry_W.insert(0, "2000")

tk.Label(frame_controls, text="Высота области H, мм").grid(row=1, column=0, sticky="w")
entry_H = tk.Entry(frame_controls, width=10)
entry_H.grid(row=1, column=1)
entry_H.insert(0, "2000")

tk.Label(frame_controls, text="Разрешение Wres, пикс").grid(row=2, column=0, sticky="w")
entry_Wres = tk.Entry(frame_controls, width=10)
entry_Wres.grid(row=2, column=1)
entry_Wres.insert(0, "400")

tk.Label(frame_controls, text="Разрешение Hres, пикс").grid(row=3, column=0, sticky="w")
entry_Hres = tk.Entry(frame_controls, width=10)
entry_Hres.grid(row=3, column=1)
entry_Hres.insert(0, "400")

# источник
tk.Label(frame_controls, text="xL, мм").grid(row=4, column=0, sticky="w")
entry_xL = tk.Entry(frame_controls, width=10)
entry_xL.grid(row=4, column=1)
entry_xL.insert(0, "0")

tk.Label(frame_controls, text="yL, мм").grid(row=5, column=0, sticky="w")
entry_yL = tk.Entry(frame_controls, width=10)
entry_yL.grid(row=5, column=1)
entry_yL.insert(0, "0")

tk.Label(frame_controls, text="zL, мм (>0)").grid(row=6, column=0, sticky="w")
entry_zL = tk.Entry(frame_controls, width=10)
entry_zL.grid(row=6, column=1)
entry_zL.insert(0, "1000")

tk.Label(frame_controls, text="I0, Вт/ср").grid(row=7, column=0, sticky="w")
entry_I0 = tk.Entry(frame_controls, width=10)
entry_I0.grid(row=7, column=1)
entry_I0.insert(0, "1000")

# круг
tk.Label(frame_controls, text="cx, мм (центр круга)").grid(row=8, column=0, sticky="w")
entry_cx = tk.Entry(frame_controls, width=10)
entry_cx.grid(row=8, column=1)
entry_cx.insert(0, "0")

tk.Label(frame_controls, text="cy, мм (центр круга)").grid(row=9, column=0, sticky="w")
entry_cy = tk.Entry(frame_controls, width=10)
entry_cy.grid(row=9, column=1)
entry_cy.insert(0, "0")

tk.Label(frame_controls, text="R, мм (радиус круга)").grid(row=10, column=0, sticky="w")
entry_R = tk.Entry(frame_controls, width=10)
entry_R.grid(row=10, column=1)
entry_R.insert(0, "800")


# Связываем поля W/H и Wres/Hres для автоматического поддержания квадратных пикселей
W_var = tk.StringVar(value="2000")
H_var = tk.StringVar(value="2000")
Wres_var = tk.StringVar(value="400")
Hres_var = tk.StringVar(value="400")

entry_W.config(textvariable=W_var)
entry_H.config(textvariable=H_var)
entry_Wres.config(textvariable=Wres_var)
entry_Hres.config(textvariable=Hres_var)

def on_W_change(*args):
    try:
        W = float(W_var.get())
        H = float(H_var.get())
        Hres = int(Hres_var.get())
        Wres_new = int(round(Hres * W / H))
        Wres_var.set(str(Wres_new))
    except:
        pass

def on_H_change(*args):
    try:
        W = float(W_var.get())
        H = float(H_var.get())
        Wres = int(Wres_var.get())
        Hres_new = int(round(Wres * H / W))
        Hres_var.set(str(Hres_new))
    except:
        pass

def on_Wres_change(*args):
    try:
        W = float(W_var.get())
        H = float(H_var.get())
        Wres = int(Wres_var.get())
        Hres_new = int(round(Wres * H / W))
        Hres_var.set(str(Hres_new))
    except:
        pass

def on_Hres_change(*args):
    try:
        W = float(W_var.get())
        H = float(H_var.get())
        Hres = int(Hres_var.get())
        Wres_new = int(round(Hres * W / H))
        Wres_var.set(str(Wres_new))
    except:
        pass

W_var.trace_add("write", on_W_change)
H_var.trace_add("write", on_H_change)
Wres_var.trace_add("write", on_Wres_change)
Hres_var.trace_add("write", on_Hres_change)


# Кнопка сохранения
def on_save():
    global last_result
    if last_result is None:
        messagebox.showwarning("Внимание", "Сначала выполните расчёт.")
        return

    E_norm = last_result["E_norm"]
    file_path = filedialog.asksaveasfilename(
        defaultextension=".png",
        filetypes=[("PNG files", "*.png"), ("All files", "*.*")]
    )
    if not file_path:
        return

    from matplotlib import pyplot as plt

    fig_save, ax = plt.subplots()
    ax.imshow(E_norm, cmap="gray", origin="lower")
    ax.axis("off")
    fig_save.savefig(file_path, bbox_inches="tight", pad_inches=0)
    plt.close(fig_save)
    messagebox.showinfo("Сохранение", f"Изображение сохранено:\n{file_path}")


btn_compute = tk.Button(frame_controls, text="Рассчитать и визуализировать")
btn_compute.grid(row=11, column=0, columnspan=2, pady=(10, 5), sticky="we")

btn_save = tk.Button(frame_controls, text="Сохранить изображение", command=on_save)
btn_save.grid(row=12, column=0, columnspan=2, pady=(0, 10), sticky="we")

# Поле для текстовых результатов
tk.Label(frame_controls, text="Результаты:").grid(row=13, column=0, columnspan=2, sticky="w")
results_text = tk.Text(frame_controls, width=40, height=12)
results_text.grid(row=14, column=0, columnspan=2, pady=(2, 0))


# ----- Фигура matplotlib внутри Tk -----

fig = Figure(figsize=(12, 9), dpi=120)
gs = GridSpec(2, 1, figure=fig, height_ratios=[1, 1])

ax_img = fig.add_subplot(gs[0])   # верхний — распределение
ax_prof = fig.add_subplot(gs[1])  # нижний — сечение

canvas = FigureCanvasTkAgg(fig, master=frame_plot)
canvas_widget = canvas.get_tk_widget()
canvas_widget.pack(fill="both", expand=True)

cbar = None
last_result = None


def update_plots(res, W_mm, H_mm, cx_mm, cy_mm, R_mm):
    global cbar

    E = res["E"]                # Вт/м^2
    E_norm = res["E_norm"]
    xs_mm = res["xs_mm"]
    ys_mm = res["ys_mm"]
    circle_mask = res["circle_mask"]

    # ----- верхний график: карта освещенности (только внутри круга видно) -----
    ax_img.clear()

    extent = [-W_mm / 2.0, W_mm / 2.0, -H_mm / 2.0, H_mm / 2.0]

    im = ax_img.imshow(
        E_norm,
        cmap="gray",
        origin="lower",
        extent=extent,
        aspect="equal"
    )

    # круг границы области визуализации
    circle = Circle((cx_mm, cy_mm), R_mm, fill=False, linewidth=1.5, color='black')
    ax_img.add_patch(circle)

    ax_img.set_title("Нормированное распределение освещенности (только внутри круга)")
    ax_img.set_xlabel("X, мм")
    ax_img.set_ylabel("Y, мм")

    # цветовая шкала
    if cbar is None:
        cbar = fig.colorbar(im, ax=ax_img, fraction=0.046, pad=0.04)
        cbar.set_label("Уровень (0-255)")
    else:
        cbar.update_normal(im)

    # ----- нижний график: горизонтальное сечение через центр круга -----
    ax_prof.clear()

    # выбираем строку, ближайшую к y = cy_mm
    center_row = np.argmin(np.abs(ys_mm - cy_mm))
    section_E = E[center_row, :]             # Вт/м^2
    section_mask = circle_mask[center_row, :]

    # только внутри круга, остальное не выводим
    section = np.where(section_mask, section_E, np.nan)

    ax_prof.plot(xs_mm, section)
    ax_prof.set_xlabel("X, мм")
    ax_prof.set_ylabel("E, Вт/м²")
    ax_prof.set_title("Освещенность вдоль горизонтального сечения (только внутри круга)")
    ax_prof.grid(True)

    canvas.draw()


def on_compute():
    try:
        W_mm = float(entry_W.get())
        H_mm = float(entry_H.get())
        Wres = int(entry_Wres.get())
        Hres = int(entry_Hres.get())
        xL = float(entry_xL.get())
        yL = float(entry_yL.get())
        zL = float(entry_zL.get())
        I0 = float(entry_I0.get())
        cx = float(entry_cx.get())
        cy = float(entry_cy.get())
        R = float(entry_R.get())
    except ValueError:
        messagebox.showerror("Ошибка", "Проверьте ввод числовых параметров.")
        return

    if zL <= 0:
        messagebox.showerror("Ошибка", "zL должен быть > 0 (источник над плоскостью).")
        return

    res = compute_illuminance(
        W_mm, H_mm,
        Wres, Hres,
        xL, yL, zL,
        I0,
        cx, cy, R
    )

    global last_result
    last_result = res

    # Обновляем графики
    update_plots(res, W_mm, H_mm, cx, cy, R)

    # Текст результатов (Вт/м^2)
    circle_stats = res["circle_stats"]
    cp = res["control_points"]

    lines = []
    lines.append("Пять контрольных точек (E, Вт/м²):")
    for name, val in cp.items():
        lines.append(f"{name}: {val:.6e}")
    lines.append("")
    lines.append("Статистика по кругу (E, Вт/м²):")
    lines.append(f"E_min: {circle_stats['E_min']:.6e}")
    lines.append(f"E_max: {circle_stats['E_max']:.6e}")
    lines.append(f"E_avg: {circle_stats['E_avg']:.6e}")

    txt = "\n".join(lines)
    print(txt)

    results_text.delete("1.0", tk.END)
    results_text.insert(tk.END, txt)


btn_compute.config(command=on_compute)

if __name__ == "__main__":
    root.mainloop()
